from pathlib import Path, PosixPath

from ase import Atoms
from ase.build import surface
from ase.io import read, write

from pseudopotential_calculator.calculations.slab import (
    SlabOptimizationParams,
    get_bulk_optimization_calculator,
)
from pseudopotential_calculator.castep import Castep, CastepConfig
from pseudopotential_calculator.hpc import (
    copy_files_to_hpc,
    prepare_calculator_with_submit_script,
    prepare_submit_all_script,
)
from pseudopotential_calculator.util import prepare_clean_directory


def _prepare_relaxation(surface: Atoms) -> None:
    data_path = Path("data/copper/slab")
    calculators = list[Castep]()
    prepare_clean_directory(data_path)

    for cutoff_energy in range(300, 750, 50):
        config = CastepConfig(data_path / f"surface_{cutoff_energy}", "surface")
        params = SlabOptimizationParams(n_k_points=11, cut_off_energy=cutoff_energy)
        calculator = get_bulk_optimization_calculator(atom, params, config)
        prepare_calculator_with_submit_script(calculator)

        calculators.append(calculator)

    prepare_submit_all_script(calculators, data_path)
    copy_files_to_hpc(data_path, PosixPath("data/copper/bulk/cutoff_energy"))


if __name__ == "__main__":
    copper_bulk = read(
        "scripts/copper/Cu.cell",
        format="castep-cell",
    )  # TODO this is to be substituted with correct parameters
    copper_slab = surface(copper_bulk, (1, 1, 0), layers=5, vacuum=40)

    write("scripts/copper/cu_slab.cell", copper_slab, format="castep-cell")

    _prepare_relaxation(copper_slab)
