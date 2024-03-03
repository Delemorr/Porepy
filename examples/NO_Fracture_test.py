import porepy as pp
import numpy as np
from porepy.applications.md_grids.domains import nd_cube_domain
from porepy.models.fluid_mass_balance import BoundaryConditionsSinglePhaseFlow
from src.equivalent_permeability import ekv

from porepy.models.constitutive_laws import SecondOrderTensorUtils

number = pp.number

k1 = 100 * 10 ** (-15)
k2 = 500 * 10 ** (-15)
P1 = 10e5
P2 = 5e5
m = 0.25
"L - Длина области"
L = 10
"z-размер ячейки"
z = 0.2
"ro - плотность"
ro = 1000
"mu - вязкость"
mu = 0.01

# ---------------------------------------------------------------------------------------------
"""Класс для изменения геометрии"""
class ModifiedGeometry:

    def set_domain(self) -> None:
        """ Определение двумерной квадратной области"""
        size = self.solid.convert_units(L, "m")
        self._domain = nd_cube_domain(2, size)

    def grid_type(self) -> str:
        """ Выбор типа сетки для нашей расчетной области."""
        return self.params.get("grid_type", "simplex")

    def meshing_arguments(self) -> dict:
        """Определение размера ячейки"""
        cell_size = self.solid.convert_units(z, "m")
        mesh_args: dict[str, float] = {"cell_size": cell_size}
        return mesh_args

# ---------------------------------------------------------------------------------------------
"Класс для изменения ГУ"
class ModifiedBC(BoundaryConditionsSinglePhaseFlow):

    def bc_type_darcy_flux(self, sd: pp.Grid) -> pp.BoundaryCondition:
        """ Назначение Дирихле западной и восточной границам. Остальные по умолчанию — Нейман."""
        bounds = self.domain_boundary_sides(sd)
        bc = pp.BoundaryCondition(sd,  bounds.west + bounds.east, "dir")
        return bc

    def bc_values_pressure(self, boundary_grid: pp.BoundaryGrid) -> np.ndarray:
        """ Назначение давления на границах."""
        bounds = self.domain_boundary_sides(boundary_grid)
        values = np.zeros(boundary_grid.num_cells)
        values[bounds.west] = self.fluid.convert_units(P1, "Pa")
        values[bounds.east] = self.fluid.convert_units(P2, "Pa")
        return values

# ---------------------------------------------------------------------------------------------
"""Изменение проницаемости"""
class Permeability(SecondOrderTensorUtils):
    def operator_to_SecondOrderTensor(self, sd: pp.Grid, operator: pp.ad.Operator, fallback_value: number,) -> pp.SecondOrderTensor:
        "Зональная неоднородность"
        x = sd.cell_centers[0]

        "Слоистая неоднородность"
        # x = sd.cell_centers[1]

        indx1 = x > L/2
        indx2 = x <= L/2
        k = np.zeros(len(x))
        k[indx1] = k1
        k[indx2] = k2
        pp.SecondOrderTensor(k, k)
        return pp.SecondOrderTensor(k, k)
# ---------------------------------------------------------------------------------------------
"""Объединение модели"""
class Flow_Model(
    Permeability,
    ModifiedGeometry,
    ModifiedBC,
    pp.fluid_mass_balance.SinglePhaseFlow,
):
    ...

# ---------------------------------------------------------------------------------------------
"""изменение параметров"""
fluid_constants = pp.FluidConstants({"viscosity": mu, "density": ro})
solid_constants = pp.SolidConstants({"porosity": m})
material_constants = {"fluid": fluid_constants, "solid": solid_constants}
params = {"material_constants": material_constants}
model = Flow_Model(params)

# ---------------------------------------------------------------------------------------------
"""Запуск модели"""
pp.run_time_dependent_model(model, params)
"Построение графика"
pp.plot_grid(model.mdg, "pressure", figsize=(10, 8), linewidth=0.25, title="Pressure distribution", plot_2d=True)

# ---------------------------------------------------------------------------------------------
"""Информация про проницаемость"""
sd = model.mdg.subdomains()[0]
data = model.mdg.subdomain_data(sd)
gg = data['parameters']['flow']['second_order_tensor'].values
pp.plot_grid(sd, gg[0,0], title="Проницаемость", plot_2d=True)

# ---------------------------------------------------------------------------------------------
"""Поиск эквивалетной проницаемости"""
k = ekv(z, L, ro, mu, P1, P2, model)
print("Эквивалетная проницаемость \n" + str(k))
