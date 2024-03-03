import numpy as np
import porepy as pp
from porepy.applications.md_grids.domains import nd_cube_domain
from src.equivalent_permeability import ekv
from porepy.models.constitutive_laws import SecondOrderTensorUtils
from porepy.applications.discretizations.flux_discretization import FluxDiscretization
from porepy.models.constitutive_laws import DimensionDependentPermeability
from typing import Callable, Union

from src.romb_fracture import romb_fractur
from src.hexagon_fracture import hexagon_fractur
from src.square_fracture import square_fractur


number = pp.number

k1 = 100 * 10 ** (-15)
k2 = 500 * 10 ** (-15)
# k1 = 1
# k2 = 1
ktr = 10000 * 10 ** (-15)
P1 = 10e5
P2 = 5e5
m = 0.25
"L - Длина области"
L = 10
"b - размер стороны трещины"
b = 2
"z-размер ячейки"
z = 0.5
"ro - плотность"
ro = 1000
"mu - вязкость"
mu = 0.01

# ---------------------------------------------------------------------------------------------
"Построение системы трещин"
box = {'xmin': 0, 'xmax': L, 'ymin': 0, 'ymax': L}
domain = pp.Domain(bounding_box=box)
fractures = romb_fractur(L,b)
network_2d = pp.create_fracture_network(fractures, domain)
mesh_args: dict[str, float] = {"cell_size": 1.0, "cell_size_fracture": 0.9}
mdg = pp.create_mdg("simplex", mesh_args, network_2d)
"Построение графика"
network_2d.plot()
# ---------------------------------------------------------------------------------------------
"Однако класс SolidConstants по умолчанию не охватывает отдельные трещины и проницаемость матрицы"

class FractureSolidConstants(pp.SolidConstants):
    """Твердые константы, адаптированные к текущей модели"""

    @property
    def default_constants(self):
        """Добавьте дополнительный параметр «fracture_permeability»."""
        constants = super().default_constants
        constants.update({"fracture_permeability": 1})
        return constants

    def fracture_permeability(self) -> float:
        """Проницаемость трещин [м^2]"""
        return self.convert_units(self.constants["fracture_permeability"], "m^2")

# ---------------------------------------------------------------------------------------------
"""Класс для изменения геометрии"""

class Geometry:

    def set_domain(self) -> None:
        """ Определение двумерной квадратной области"""
        size = self.solid.convert_units(L, "m")
        self._domain = nd_cube_domain(2, size)

    def set_fractures(self) -> None:
        """Установка Трещин"""
        self._fractures = romb_fractur(L, b)

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
class BoundaryConditions:
    domain_boundary_sides: Callable[[pp.Grid | pp.BoundaryGrid], pp.domain.DomainSides]
    """Граничные стороны области. Определяется экземпляром примеси
    :class:`~porepy.models.geometry.ModelGeometry`.

    """
    fluid: pp.FluidConstants
    """Объект константы жидкости, который обеспечивает масштабирование величин, связанных с жидкостью.
     Обычно это устанавливается миксином экземпляра
    :class:`~porepy.models.solution_strategy.SolutionStrategy`.

    """
    specific_volume: Callable[
        [Union[list[pp.Grid], list[pp.MortarGrid]]], pp.ad.Operator]
    """Функция, возвращающая конкретный объем поддомена или интерфейса.
     Обычно предоставляется примесью экземпляра
    :class:`~porepy.models.constitutive_laws.DimensionReduction`.

    """
    equation_system: pp.ad.EquationSystem
    """Объект EquationSystem для текущей модели. Обычно определяется в классе примеси
     определение стратегии решения.

    """
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

    def bc_values_darcy_flux(self, boundary_grid: pp.BoundaryGrid) -> np.ndarray:
        """Приток на западной границе.

         Согласно соглашению PorePy, знак притока отрицательный, а значение равно
         интегрированный по объемам граничных ячеек. Поскольку граница притока содержит
         трещины, последний включает удельный объем трещины.

         Параметры:
             border_grid: Граничная сетка.

         Возврат:
             Граничные значения.

        """
        bounds = self.domain_boundary_sides(boundary_grid)
        values = np.zeros(boundary_grid.num_cells)
        # Приток на западной границе. Подпишите согласно соглашению PorePy
        val = self.fluid.convert_units(-1, "m * s^-1")
        # Интегрируйте по объемам граничных ячеек.
        values[bounds.west] = val * boundary_grid.cell_volumes[bounds.west]
        # Масштабируйте с определенным объемом.
        sd = boundary_grid.parent
        trace = np.abs(sd.cell_faces)
        specific_volumes = self.specific_volume([sd]).value(self.equation_system)
        values *= boundary_grid.projection() @ trace @ specific_volumes
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

class Permeabilit(DimensionDependentPermeability):
    params: dict
    solid: FractureSolidConstants

    def fracture_permeability(self, subdomains: list[pp.Grid]) -> pp.ad.Operator:
        size = sum([sd.num_cells for sd in subdomains])
        permeability = pp.wrap_as_dense_ad_array(
            self.solid.fracture_permeability(), size, name="fracture permeability"
        )
        return self.isotropic_second_order_tensor(subdomains, permeability)

# ---------------------------------------------------------------------------------------------
"""Объединение модели"""
class Flow_Model(
    FluxDiscretization,
    Permeabilit,
    Geometry,
    Permeability,
    BoundaryConditions,
    pp.fluid_mass_balance.SinglePhaseFlow,):
    ...
# ---------------------------------------------------------------------------------------------
"""изменение параметров для системы трещин"""

solid_constants_conductive_fractures = FractureSolidConstants(
    {
        "residual_aperture": 1e-2,
        "fracture_permeability": ktr,
        "normal_permeability": ktr,
    }
)
# solid_constants_blocking_fractures = FractureSolidConstants(
#     {
#         # "residual_aperture": 1e-4,
#         "fracture_permeability": 1,
#         "normal_permeability": 1,
#     }
# )

"Модуль содержит реализацию схемы аппроксимации двухточечного потока конечного объема. Реализация находится в классе Tpfa."

# solid_constants = [solid_constants_blocking_fractures, solid_constants_conductive_fractures]
# solid_constants = [solid_constants_blocking_fractures]
solid_constants = [solid_constants_conductive_fractures]
for solid_constants, discr in zip(solid_constants, ["tpfa", "mpfa"]):
    # Мы используем параметры жидкости по умолчанию, но адаптированные параметры твердого тела.
    fluid_constants = pp.FluidConstants({"viscosity": mu, "density": ro})
    # solid_constant = pp.SolidConstants({"porosity": m})
    params = {"material_constants": {"fluid": fluid_constants, "solid": solid_constants}}

# ---------------------------------------------------------------------------------------------
"""Запуск модели"""
model = Flow_Model(params)
pp.run_time_dependent_model(model, params)
pp.plot_grid(model.mdg, "pressure", figsize=(10, 8), linewidth=0.25, title="Pressure distribution", plot_2d=True, pointsize=20, fracturewidth_1d=3)

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


