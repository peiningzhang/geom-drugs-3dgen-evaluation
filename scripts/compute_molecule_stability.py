import argparse

import numpy as np
from rdkit import Chem
from rdkit import RDLogger

# Disable all RDKit warnings
RDLogger.DisableLog('rdApp.*')

from geom_utils.molecule_stability import compute_molecules_stability


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--sdf", type=str, required=True, help="Path to input SDF file")
    parser.add_argument("--n_subsets", type=int, default=1,
                        help="Number of subsets to compute std over.")
    parser.add_argument("--aromatic", type=str, choices=["true", "false", "auto"], default="auto",
                        help="Whether to consider aromatic bonds: true, false, or auto")

    args = parser.parse_args()

    mols = [m for m in Chem.SDMolSupplier(args.sdf, sanitize=False, removeHs=False) if
            m is not None]

    if args.aromatic == "true":
        aromatic = True
    elif args.aromatic == "false":
        aromatic = False
    else:
        aromatic = any(b.GetIsAromatic() for m in mols for b in m.GetBonds())

    if args.n_subsets is None:
        validity, stable_mask, n_stable_atoms, n_atoms = compute_molecules_stability(mols,
                                                                                     aromatic=aromatic)
        print(f"Valid molecules: {int(validity.sum())} / {len(mols)} ({validity.mean():.4f})")
        print(
            f"Stable molecules: {int(stable_mask.sum())} / {len(mols)} ({stable_mask.mean():.4f})")
        print(f"Total atoms: {int(n_atoms.sum())}")
        print(
            f"Stable atoms: {int(n_stable_atoms.sum())} ({(n_stable_atoms.sum() / n_atoms.sum()):.4f})")
    elif args.n_subsets == -1:
        # 统计各个atom_number下同时Valid和stable的比例
        from collections import defaultdict
        validity, stable_mask, n_stable_atoms, n_atoms = compute_molecules_stability(mols,
                                                                                     aromatic=aromatic)
        atom_count_stats = defaultdict(lambda: {'total': 0, 'valid_and_stable': 0})
        
        for i, mol in enumerate(mols):
            if mol is not None:
                atom_count = n_atoms[i].item()
                is_valid = validity[i].item() > 0
                is_stable = stable_mask[i].item() > 0
                
                atom_count_stats[atom_count]['total'] += 1
                if is_valid and is_stable:
                    atom_count_stats[atom_count]['valid_and_stable'] += 1
        
        # 按atom_number从小到大排序并打印
        print("\nMolecules that are both Valid and Stable by atom count:")
        for atom_count in sorted(atom_count_stats.keys()):
            stats = atom_count_stats[atom_count]
            ratio = stats['valid_and_stable'] / stats['total'] if stats['total'] > 0 else 0
            print(f"Atom count {atom_count}: {stats['valid_and_stable']} / {stats['total']} ({ratio:.4f})")
        print(f"Valid molecules: {int(validity.sum())} / {len(mols)} ({validity.mean():.4f})")
        print(
            f"Stable molecules: {int(stable_mask.sum())} / {len(mols)} ({stable_mask.mean():.4f})")
        print(f"Total atoms: {int(n_atoms.sum())}")
        print(
            f"Stable atoms: {int(n_stable_atoms.sum())} ({(n_stable_atoms.sum() / n_atoms.sum()):.4f})")
    else:
        mols_per_set = len(mols) // args.n_subsets
        valids, mols_stable, atoms_stable = [], [], []

        for i in range(args.n_subsets):
            chunk = mols[i * mols_per_set:(i + 1) * mols_per_set]
            validity, stable_mask, n_stable_atoms, n_atoms = compute_molecules_stability(chunk,
                                                                                         aromatic=aromatic)
            valids.append(validity.mean().item())
            mols_stable.append(stable_mask.mean().item())
            atoms_stable.append((n_stable_atoms.sum() / n_atoms.sum()).item())

        print(f"Avg validity: {np.mean(valids):.4f} ± {np.std(valids):.4f}")
        print(f"Avg molecule stability: {np.mean(mols_stable):.4f} ± {np.std(mols_stable):.4f}")
        print(f"Avg atom stability: {np.mean(atoms_stable):.4f} ± {np.std(atoms_stable):.4f}")


if __name__ == "__main__":
    main()
