#!/usr/bin/env python3
"""
Recombination Benchmarker

Benchmarks the native GRG recombination approach against the baseline
NumPy conversion approach (Recombine -> Dense Matrix).
"""

import os
import glob
import argparse
import time
import csv
import json
import gc
import sys
import platform
import numpy as np
from pathlib import Path
from dataclasses import dataclass, asdict
import pygrgl
import re
import psutil
import cProfile, pstats

# Import our custom modules (ensure these are in the same directory)
from grg_recombination import simulate_grg_recombination, NonDuplicationRecombination
from grg_numpy_baseline import grg_to_numpy, grg_to_numpy_parallel, estimate_numpy_memory, simulate_numpy_recombination

# DEBUG Mode
debug = True

def get_process_memory_mb():
    """Returns the current process memory usage in MB."""
    process = psutil.Process(os.getpid())
    # rss = Resident Set Size (the actual RAM used)
    return process.memory_info().rss / (1024 * 1024)

# def get_breakpoints(bp_range, expected_crossovers=1.5):
#     num_bp = np.random.poisson(expected_crossovers)
#     if num_bp == 0:
#         return np.array([], dtype=int)

#     low, high = bp_range            # inclusive or exclusive? assuming inclusive high here
#     length = high - low + 1
#     num_bp = min(num_bp, length)

#     # Oversample a bit, then deduplicate
#     k = num_bp * 3                  # oversampling factor; small, since num_bp is tiny
#     candidates = np.random.randint(low, high + 1, size=k)
#     unique = np.unique(candidates)

#     if unique.size < num_bp:
#         # Very unlikely with small num_bp, but handle just in case
#         extra_needed = num_bp - unique.size
#         extra = np.random.randint(low, high + 1, size=extra_needed * 3)
#         unique = np.unique(np.concatenate([unique, extra]))

#     bp = np.sort(unique[:num_bp])
#     return bp


def recombination_intervals(h1, h2, bp, N):
    """
    Returns list of segments: [(source_parent_id, end_coord), ...]
    """
    # bp = self.get_breakpoints(N)
    # start = np.random.binomial(1, 0.5, 1)[0]
    start = h1
    parents = [h1, h2]
    
    segments = []
    for i, K in enumerate(bp):
        segments.append((parents[(start + i) % 2], K))
    segments.append((parents[(start + len(bp)) % 2], N))
    return segments

@dataclass
class BenchmarkResult:
    """Single benchmark measurement result"""
    file: str
    num_offspring: int
    implementation: str  # 'GRG Native' or 'NumPy Baseline'
    num_samples_initial: int
    num_snps: int
    num_runs: int
    num_bp: float
    mean_bp: float
    mean_time_ms: float
    std_time_ms: float
    min_time_ms: float
    max_time_ms: float
    nodes_added: int = 0
    memory_mb: float = 0.0

@dataclass
class SystemInfo:
    """System and environment information"""
    numpy_version: str
    pygrgl_version: str
    platform: str
    cpu_count: int
    timestamp: str

class RecombinationBenchmarker:
    def __init__(
        self,
        grg_files_dir: Path,
        output_dir: Path,
        num_warmup: int = 1,
        num_runs: int = 3,
        num_offspring_per_couple: int = 2,
        num_generations: int = 2,
        memory_limit_mb: float = 15000.0,
        include_numpy: bool = True
    ):
        self.grg_dir = grg_files_dir
        self.output_dir = output_dir
        self.num_warmup = num_warmup
        self.num_runs = num_runs
        self.num_offspring_per_couple = num_offspring_per_couple
        self.num_couples = None  # will be determined per file based on num_samples
        self.num_generations = num_generations
        self.memory_limit_mb = memory_limit_mb
        self.include_numpy = include_numpy
        self.results = []
        
        try:
            pygrgl_ver = pygrgl.__version__
        except AttributeError:
            pygrgl_ver = "unknown"
            
        import multiprocessing
        self.system_info = SystemInfo(
            numpy_version=np.__version__,
            pygrgl_version=pygrgl_ver,
            platform=f"{platform.system()}-{platform.release()}-{platform.machine()}",
            cpu_count=multiprocessing.cpu_count(),
            timestamp=time.strftime("%Y-%m-%d %H:%M:%S")
        )

    def _parse_expected_size(self, filename: str) -> tuple:
        """Extract expected individuals and SNPs from filename."""
        inds_match = re.search(r'(\d+)inds', filename)
        snps_match = re.search(r'(\d+)([kKmM])?_snps', filename)
        
        inds = int(inds_match.group(1)) if inds_match else 0
        
        snps = 0
        if snps_match:
            val = int(snps_match.group(1))
            mult = snps_match.group(2)
            if mult:
                mult = mult.lower()
                if mult == 'k': snps = val * 1000
                elif mult == 'm': snps = val * 1000000
            else:
                snps = val
                
        return inds, snps

    def print_header(self):
        print("=" * 80)
        print("GRG Recombination vs NumPy Baseline Benchmark")
        print("=" * 80)
        print(f"System: {self.system_info.platform}")
        print(f"CPU Cores: {self.system_info.cpu_count}")
        print(f"NumPy: {self.system_info.numpy_version}")
        print(f"Warmup runs: {self.num_warmup}, Timed runs: {self.num_runs}")
        print(f"Offspring per couple: {self.num_offspring_per_couple}")
        print(f"Memory limit: {self.memory_limit_mb} MB")
        print("=" * 80)

    def get_breakpoints(self, bp_range, expected_crossovers=1.5):
        num_bp = np.random.poisson(expected_crossovers)
        if num_bp == 0:
            return np.array([], dtype=int), num_bp

        low, high = bp_range  # assuming exclusive low and high
        length = high - low + 1
        num_bp = min(num_bp, length)

        # Oversample a bit, then deduplicate
        k = num_bp * 3  # oversampling factor; small, since num_bp is tiny
        candidates = np.random.randint(low + 1, high, size=k)
        unique = np.unique(candidates)

        if unique.size < num_bp:
            # Very unlikely with small num_bp, but handle just in case
            extra_needed = num_bp - unique.size
            extra = np.random.randint(low + 1, high, size=extra_needed * 3)
            unique = np.unique(np.concatenate([unique, extra]))

        bp = np.sort(unique[:num_bp])
        return bp, num_bp

    def run_benchmarks(self):
        self.print_header()

        # 1. Use pathlib to find the files (returns real Path objects, not strings!)
        if self.grg_dir.is_dir():
            grg_files = list(self.grg_dir.glob("*.grg"))
        else:
            grg_files = [self.grg_dir] if self.grg_dir.suffix == ".grg" else []

        # 2. Sort from smallest to largest using Path's native stat().st_size
        grg_files.sort(key=lambda x: x.stat().st_size)
        
        if not grg_files:
            print(f"No .grg files found in {self.grg_dir}")
            return []

        print("Files to process (Smallest to Largest):")
        for f in grg_files:
            # f is now a real Path object, so .name works perfectly here!
            print(f"  {f.name}")
        print("=" * 80)

        # Set base csv output filename prefix (we'll append file-specific suffix if benchmark is run on 1 file)
        filename_prefix = "benchmark_recombination_results"
        if len(grg_files) == 1:
            # If only one file, include its name in the output filename prefix for clarity
            filename_prefix += f"_{grg_files[0].stem}"

        for file_path in grg_files:
            # file_path is a real Path object, so _benchmark_file will accept it perfectly!
            self._benchmark_file(file_path)
            
            self.save_results(self.output_dir, filename_prefix)
            print(f"[*] Checkpoint saved! Results safely written to disk.")
        return self.results

    def _benchmark_file(self, file_path: Path):
        print(f"\nProcessing: {file_path.name}")
        exp_inds, exp_snps = self._parse_expected_size(file_path.name)
        if exp_inds and exp_snps:
            print(f" Expected: {exp_inds} individuals, {exp_snps} SNPs")
            
        # Load once just to get base stats
        g_base = pygrgl.load_mutable_grg(str(file_path))
        base_samples = g_base.num_samples
        base_mutations = g_base.num_mutations
        base_nodes = g_base.num_nodes
        base_genome = g_base.bp_range
        base_couples = base_samples // 2
        self.num_couples = base_couples
        
        print(f" Genome (bp): {base_genome}")
        print(f" Couples: {base_couples} (from {base_samples} samples)")

        # Estimate memory for the POST-recombination NumPy array
        # It will have (base_samples + num_offspring) columns
        est_mb = estimate_numpy_memory(base_mutations, base_samples + self.num_offspring_per_couple * base_couples)

        skip_numpy = False
        if est_mb > self.memory_limit_mb:
            print(f"  [Skipping NumPy] Estimated memory ({est_mb:.1f} MB) exceeds limit ({self.memory_limit_mb} MB)")
            skip_numpy = True

        print(f" Actual Base: {base_samples} samples, {base_mutations} mutations, {base_nodes} nodes")
        print(f" Estimated memory for NumPy array (post-recomb): {est_mb:.1f} MB")
        
        # Free base graph
        del g_base
        gc.collect()

        # generate segments
        # run grg recombination
        # run numpy recombination
        # repeat steps 1-3 for num_runs and num_warmup

        print(f"\nGenerating Base NumPy Population Matrix...")
        base_numpy_population_matrix = grg_to_numpy_parallel(pygrgl.load_mutable_grg(str(file_path)))

        grg_times = []
        grg_sizes = []
        grg_size_changes = []
        numpy_times = []
        nodes_added = 0

        for i in range(self.num_warmup + self.num_runs):

            print(f"  [Run {i+1}] Loading graph from disk...")
            load_start = time.time()
            g = pygrgl.load_mutable_grg(str(file_path))
            print(f"  [Run {i+1}] Graph loaded in {time.time() - load_start:.2f} seconds.")

            # testcases = []

            # # Generate testcase breakpoints and segments for all offspring in this run (for both GRG and NumPy)
            # if debug:
            #     print(f"  [Run {i+1}] Generating testcases")
            # for gen in range(self.num_generations):
            #     gen_testcases = []
            #     for couple in range(base_couples):
            #         couple_testcases = []
            #         for offspring in range(4):  # Each pair produces 2 children and each child has 2 sets of breakpoints due to diploid instance
            #             print("testcase breakpoints:", couple, offspring)
            #             bp = get_breakpoints(base_genome)
            #             couple_testcases.append(bp)
            #         gen_testcases.append(couple_testcases)
            #     testcases.append(gen_testcases)
            # if debug:
            #     print(f"  [Run {i+1}] Testcases generated.")



            # GRG Native
            gc.collect()
            mem_before = get_process_memory_mb()
            print(f"\nStarting GRG Recombination Benchmarking...")
            if debug:
                print(f"\n  [Run {i+1}] Simulating {self.num_generations} generations with GRG recombination...")
            total_grg_bp = 0
            start = time.perf_counter()
            recomb = NonDuplicationRecombination(g)
            for gen in range(self.num_generations):
                # Run GRG recombination for this generation
                if debug:
                    print(f"\n    [Gen {gen+1}] Simulating generation {gen+1} with GRG recombination...")
                # pr = cProfile.Profile()
                # pr.enable()
                offspring_ids, gen_bp = simulate_grg_recombination(self, recomb, base_genome, N=base_genome[1])
                # pr.disable()
                # pstats.Stats(pr).sort_stats('cumtime').print_stats(25)  # Print top 25 cumulative time functions
                total_grg_bp += gen_bp
                if debug:
                    print(f"    [Gen {gen+1}] Generation's Breakpoints: {gen_bp}. Total so far: {total_grg_bp}")
            elapsed = time.perf_counter() - start
            if debug:
                print(f"\n  [Run {i+1}] Simulation finished in {elapsed:.4f} seconds.")
            
            if i >= self.num_warmup:
                grg_times.append(elapsed * 1000)  # Convert to ms

            # Grab node count on the last run for space complexity metrics
            if i == (self.num_warmup + self.num_runs - 1):
                nodes_added = g.num_nodes - base_nodes
            
            gc.collect()
            mem_after = get_process_memory_mb()
            grg_sizes.append(mem_after)
            grg_size_changes.append(mem_after - mem_before)
            print(f"  [Run {i+1}] Memory usage: {mem_after:.1f} MB (Delta: +{mem_after - mem_before:.1f} MB)")
            del g
            del recomb

            gc.collect()
            # Numpy
            if not skip_numpy:
                print(f"\nStarting NumPy Recombination Benchmarking...")
                if debug:
                    print(f"\n  [Run {i+1}] Simulating {self.num_generations} generations with NumPy recombination...")
                offspring_matrix = base_numpy_population_matrix.copy()
                total_numpy_bp = 0
                start = time.perf_counter()
                for gen in range(self.num_generations):
                    # Run NumPy recombination for this generation
                    if debug:
                        print(f"\n    [Gen {gen+1}] Simulating generation {gen+1} with NumPy recombination...")
                    offspring_matrix, gen_bp = simulate_numpy_recombination(self, offspring_matrix, base_genome)
                    total_numpy_bp += gen_bp
                    if debug:
                        print(f"    [Gen {gen+1}] Generation's Breakpoints: {gen_bp}. Total so far: {total_numpy_bp}")
                elapsed = time.perf_counter() - start
                if debug:
                    print(f"\n  [Run {i+1}] Simulation finished in {elapsed:.4f} seconds.")
                if i >= self.num_warmup:
                    numpy_times.append(elapsed * 1000)  # Convert to ms
        
        grg_mean, grg_std, grg_sizes_mean, grg_size_changes_mean = np.mean(grg_times), np.std(grg_times), np.mean(grg_sizes), np.mean(grg_size_changes)
        grg_bp_mean = total_grg_bp / (self.num_runs * self.num_generations)
        if not skip_numpy:
            numpy_mean, numpy_std = np.mean(numpy_times), np.std(numpy_times)
            numpy_bp_mean = total_numpy_bp / (self.num_runs * self.num_generations)
        actual_offspring_generated = self.num_couples * self.num_offspring_per_couple * self.num_generations

        self.results.append(BenchmarkResult(
            file=file_path.name, num_offspring=actual_offspring_generated,
            implementation="GRG Native", num_samples_initial=base_samples,
            num_snps=base_mutations, num_runs=self.num_runs,
            num_bp=total_grg_bp, mean_bp=grg_bp_mean,
            mean_time_ms=grg_mean, std_time_ms=grg_std,
            min_time_ms=np.min(grg_times), max_time_ms=np.max(grg_times),
            nodes_added=nodes_added, memory_mb=grg_sizes_mean,
        ))

        if not skip_numpy:
            self.results.append(BenchmarkResult(
                file=file_path.name, num_offspring= actual_offspring_generated,
                implementation="NumPy Baseline", num_samples_initial=base_samples,
                num_snps=base_mutations, num_runs=self.num_runs,
                num_bp=total_numpy_bp, mean_bp=numpy_bp_mean,
                mean_time_ms=numpy_mean, std_time_ms=numpy_std,
                min_time_ms=np.min(numpy_times), max_time_ms=np.max(numpy_times),
                nodes_added=0, memory_mb=est_mb
            ))

        print(f"\nResults for {file_path.name}:\n")

        print(f"  GRG Native:     {grg_mean:.2f}ms ± {grg_std:.2f}ms")
        print(f"  Space Delta:    +{nodes_added} nodes created")
        print(f"  GRG Breakpoints: {grg_bp_mean:.2f} average breakpoints per generation")
        print(f"  GRG Memory:     {grg_sizes_mean:.1f} MB average resident memory")
        print(f"  Memory Delta:   +{grg_size_changes_mean:.1f} MB (GRG recombination memory increase), ")
        print()

        if not skip_numpy:
            print(f"  NumPy Baseline: {numpy_mean:.2f}ms ± {numpy_std:.2f}ms")
            print(f"  NumPy Breakpoints: {numpy_bp_mean:.2f} average breakpoints per generation")
            if grg_mean > 0:
                speedup = numpy_mean / grg_mean
                print(f"  Speedup:        {speedup:.2f}x (NumPy / GRG)")
            print(f"  Memory Footprint: {est_mb:.1f} MB dense matrix allocated")
        
        # # ---------------------------------------------------------
        # # Phase A: GRG Native Recombination
        # # ---------------------------------------------------------

        # print("Starting GRG Benchmarking...")

        # grg_times = []
        # nodes_added = 0
        
        # for i in range(self.num_warmup + self.num_runs):
        #     # We MUST reload the graph fresh every run so it doesn't infinitely grow
        #     print(f"  [Run {i+1}] Loading graph from disk...")
        #     load_start = time.time()
        #     g = pygrgl.load_mutable_grg(str(file_path))
        #     print(f"  [Run {i+1}] Graph loaded in {time.time() - load_start:.2f} seconds.")
            
        #     print(f"  [Run {i+1}] Simulating {self.num_offspring} offspring...")
        #     start = time.perf_counter()
        #     simulate_generation(g)
        #     elapsed = time.perf_counter() - start
        #     print(f"  [Run {i+1}] Simulation finished in {elapsed:.4f} seconds.")
            
        #     if i >= self.num_warmup:
        #         grg_times.append(elapsed * 1000)  # Convert to ms
                
        #     # Grab node count on the last run for space complexity metrics
        #     if i == (self.num_warmup + self.num_runs - 1):
        #         nodes_added = g.num_nodes - base_nodes
                
        #     del g
        #     gc.collect()

        # grg_mean, grg_std = np.mean(grg_times), np.std(grg_times)

        # actual_offspring_generated = (g_base.num_samples // 2) * 2
        
        # self.results.append(BenchmarkResult(
        #     file=file_path.name, num_offspring=actual_offspring_generated,
        #     implementation="GRG Native", num_samples_initial=base_samples,
        #     num_snps=base_mutations, num_runs=self.num_runs,
        #     mean_time_ms=grg_mean, std_time_ms=grg_std,
        #     min_time_ms=np.min(grg_times), max_time_ms=np.max(grg_times),
        #     nodes_added=nodes_added, memory_mb=0.0
        # ))
        
        # print(f"  GRG Native:     {grg_mean:.2f}ms ± {grg_std:.2f}ms")
        # print(f"  Space Delta:    +{nodes_added} nodes created")

        # # ---------------------------------------------------------
        # # Phase B: NumPy Baseline
        # # ---------------------------------------------------------

        # print("Starting NumPy Baseline Benchmarking...")

        # # 1. SETUP (Not timed) - Convert GRG to Dense Matrix once
        # # This creates the starting population matrix (num_samples x num_snps)
        # numpy_population_matrix = grg_to_numpy(g)
        
        # numpy_times = []
        
        # for run in range(self.num_runs):
        #     gc.collect()
        #     start_time = time.perf_counter()
            
        #     # 2. EXECUTION (Timed) - Perform array-based recombination
        #     offspring_matrix = simulate_generation_numpy(numpy_population_matrix)
            
        #     end_time = time.perf_counter()
        #     numpy_times.append((end_time - start_time) * 1000)

        # np_mean = np.mean(numpy_times)
        # print(f"  NumPy Baseline: {np_mean:.2f}ms ± {np.std(numpy_times):.2f}ms")
        
        # # Calculate total memory of the baseline (initial matrix + offspring matrix)
        # total_np_memory = (numpy_population_matrix.nbytes + offspring_matrix.nbytes) / (1024 * 1024)
        # print(f"  Memory Footprint: {total_np_memory:.2f} MB dense matrices allocated")

        # if not self.include_numpy:
        #     return

        # if est_mb > self.memory_limit_mb:
        #     print(f"  [Skipping NumPy] Estimated memory ({est_mb:.1f} MB) exceeds limit ({self.memory_limit_mb} MB)")
        #     return

        # np_times = []
        # for i in range(self.num_warmup + self.num_runs):
        #     g = pygrgl.load_mutable_grg(str(file_path))
            
        #     # The entire baseline process: Recombine + Materialize to Dense Matrix
        #     start = time.perf_counter()
        #     simulate_generation(g, self.num_offspring)
        #     _ = grg_to_numpy(g)
        #     elapsed = time.perf_counter() - start
            
        #     if i >= self.num_warmup:
        #         np_times.append(elapsed * 1000)
                
        #     del g
        #     gc.collect()

        # np_mean, np_std = np.mean(np_times), np.std(np_times)
        
        # self.results.append(BenchmarkResult(
        #     file=file_path.name, num_offspring=self.num_offspring,
        #     implementation="NumPy Baseline", num_samples_initial=base_samples,
        #     num_snps=base_mutations, num_runs=self.num_runs,
        #     mean_time_ms=np_mean, std_time_ms=np_std,
        #     min_time_ms=np.min(np_times), max_time_ms=np.max(np_times),
        #     nodes_added=0, memory_mb=est_mb
        # ))

        # print(f"  NumPy Baseline: {np_mean:.2f}ms ± {np_std:.2f}ms")
        # if grg_mean > 0:
        #     speedup = np_mean / grg_mean
        #     print(f"  Speedup:        {speedup:.2f}x (NumPy / GRG)")
        # print(f"  Memory Footprint: {est_mb:.1f} MB dense matrix allocated")


    def save_results(self, output_dir: Path, filename_prefix: str = "benchmark_recombination_results"):
        output_dir.mkdir(parents=True, exist_ok=True)
        
        csv_path = output_dir / f"{filename_prefix}.csv"
        with open(csv_path, 'w', newline='') as f:
            if not self.results: return
            writer = csv.DictWriter(f, fieldnames=asdict(self.results[0]).keys())
            writer.writeheader()
            for res in self.results:
                writer.writerow(asdict(res))
                
        json_path = output_dir / f"{filename_prefix}.json"
        with open(json_path, 'w') as f:
            json.dump({
                "system_info": asdict(self.system_info),
                "config": {
                    "warmup_runs": self.num_warmup,
                    "timed_runs": self.num_runs,
                    "num_offspring": self.num_offspring_per_couple * self.num_couples * self.num_generations,
                    "memory_limit_mb": self.memory_limit_mb
                },
                "results": [asdict(r) for r in self.results]
            }, f, indent=2)
            
        print(f"\nResults saved to {csv_path} and {json_path}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="GRG Recombination Benchmarker")
    parser.add_argument("--grg-dir", type=Path, default=Path("./grg_files"), help="Directory with GRG files")
    parser.add_argument("--output-dir", type=Path, default=Path("."), help="Directory to save benchmark results")
    parser.add_argument("--warmup", type=int, default=1, help="Warmup runs (default: 1)")
    parser.add_argument("--runs", type=int, default=3, help="Timed runs (default: 3)")
    parser.add_argument("--offspring-per-couple", type=int, default=2, help="Number of sequential recombinations per couple in generation (default: 2)")
    parser.add_argument("--num-generations", type=int, default=2, help="Number of sequential recombination generations to simulate (default: 2)")
    parser.add_argument("--memory-limit", type=float, default=15000.0, help="Memory limit in MB (default: 15000)")
    
    args = parser.parse_args()
    
    benchmarker = RecombinationBenchmarker(
        grg_files_dir=args.grg_dir,
        output_dir=args.output_dir,
        num_warmup=args.warmup,
        num_runs=args.runs,
        num_offspring_per_couple=args.offspring_per_couple,
        num_generations=args.num_generations,
        memory_limit_mb=args.memory_limit
    )
    
    benchmarker.run_benchmarks()