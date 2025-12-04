# Multi-Process Energy Benchmark Enhancement

This enhancement adds parallel processing capabilities to the [geom-drugs-3dgen-evaluation](https://github.com/isayevlab/geom-drugs-3dgen-evaluation) benchmark repository.

## 📋 Overview

This enhancement provides:
- **Multi-process xTB optimization** for faster processing of large molecule datasets
- **Automated pipeline script** that runs the complete energy benchmark workflow

## 📦 Installation

1. Clone the original repository:
   ```bash
   git clone https://github.com/isayevlab/geom-drugs-3dgen-evaluation.git
   cd geom-drugs-3dgen-evaluation
   ```

2. Add the enhancement files:
   - Place `run_energy_benchmark.sh` in the **root directory** of the repository
   - Place `xtb_optimization_multi_process.py` in `scripts/energy_benchmark/` directory

   ```bash
   # Copy files to the repository
   cp run_energy_benchmark.sh /path/to/geom-drugs-3dgen-evaluation/
   cp xtb_optimization_multi_process.py /path/to/geom-drugs-3dgen-evaluation/scripts/energy_benchmark/
   ```

3. Make the script executable:
   ```bash
   chmod +x run_energy_benchmark.sh
   ```

## 🚀 Usage

Run the complete energy benchmark pipeline:

```bash
bash -x run_energy_benchmark.sh path/to/your/molecules.sdf
```


If no input file is specified, it defaults to `./prediction.sdf`.

## 📝 What It Does

The `run_energy_benchmark.sh` script automatically runs:

1. **xTB Optimization** (multi-process) - Optimizes molecular geometries using GFN2-xTB
2. **RMSD & Energy Analysis** - Computes energy gains and RMSD values
3. **Geometry Comparison** - Analyzes bond lengths, angles, and torsion changes
4. **Stability Assessment** - Evaluates molecular and atomic stability

All outputs are saved to timestamped log files in the `logs/` directory.

## ⚙️ Configuration

The script uses all available CPU cores by default. To limit the number of processes, modify the `xtb_optimization_multi_process.py` call in `run_energy_benchmark.sh`:

```bash
python -u scripts/energy_benchmark/xtb_optimization_multi_process.py \
    --input_sdf "$INPUT_SDF" \
    --output_sdf optimized_output.sdf \
    --init_sdf saved_initial_structures.sdf \
    --num_processes 8 \  # Specify number of processes
    >> "$LOG_FILE" 2>&1
```

## 📄 Output Files

- `optimized_output.sdf` - Optimized molecular structures
- `saved_initial_structures.sdf` - Original input structures
- `logs/YYYYMMDD_HHMMSS_energy_pipeline.log` - Complete pipeline log

## 🔧 Requirements

Same as the original repository:
- Python packages: `rdkit`, `numpy`, `tqdm`, `torch`
- External tool: `xtb` (GFN2-xTB) - must be installed and in `PATH`
- For multi-processing: `multiprocessing` (included in Python standard library)

## 📖 Original Repository

This enhancement is based on:
- Repository: [isayevlab/geom-drugs-3dgen-evaluation](https://github.com/isayevlab/geom-drugs-3dgen-evaluation)
- Paper: "GEOM-Drugs Revisited: Toward More Chemically Accurate Benchmarks for 3D Molecule Generation"
