"""Run an illustrative two-point dual-fuel nozzle study."""

from fuelnozzle import (
    JetAPropertyTable,
    LNGComposition,
    LNGFeedLine,
    LNGNozzleGeometry,
    MissionFuelMasses,
    OperatingPoint,
    PressureSwirlGeometry,
    RelaxationFlowSettings,
    run_nozzle_study,
)


def operating_point(
    name: str,
    duration_s: float,
    p3_pa: float,
    t3_k: float,
    lng_flow_kg_s: float,
    jet_a_flow_kg_s: float,
) -> OperatingPoint:
    nozzle_drop_pa = 9.0e5
    pump_pressure_pa = p3_pa + nozzle_drop_pa + 2.0e5
    return OperatingPoint(
        name=name,
        duration_s=duration_s,
        flow_multiplier=2.0,
        p3_pa=p3_pa,
        t3_k=t3_k,
        lng_mass_flow_kg_s=lng_flow_kg_s,
        jet_a_mass_flow_kg_s=jet_a_flow_kg_s,
        lng_pump_outlet_pressure_pa=pump_pressure_pa,
        lng_pump_outlet_temperature_k=120.0,
        lng_nozzle_pressure_drop_pa=nozzle_drop_pa,
        jet_a_pump_outlet_pressure_pa=pump_pressure_pa,
        jet_a_nozzle_inlet_temperature_k=300.0,
        jet_a_nozzle_pressure_drop_pa=nozzle_drop_pa,
    )


points = [
    operating_point("takeoff", 120.0, 1.2e6, 780.0, 0.12, 0.14),
    operating_point("cruise", 3_600.0, 1.0e5, 650.0, 0.08, 0.06),
]

result = run_nozzle_study(
    points,
    LNGComposition(
        mole_fractions={"Methane": 0.90, "Ethane": 0.07, "Nitrogen": 0.03}
    ),
    LNGNozzleGeometry(
        number_of_orifices=4,
        orifice_length_m=1.0e-3,
        discharge_coefficient=0.80,
    ),
    PressureSwirlGeometry(
        number_of_inlet_ports=4,
        inlet_port_diameter_m=0.7e-3,
        inlet_tangency_radius_m=2.0e-3,
        swirl_chamber_radius_m=2.5e-3,
        swirl_chamber_length_m=5.0e-3,
    ),
    JetAPropertyTable(
        temperature_k=(280.0, 300.0, 320.0),
        density_kg_m3=(815.0, 800.0, 785.0),
        viscosity_pa_s=(2.1e-3, 1.5e-3, 1.1e-3),
        surface_tension_n_m=(0.027, 0.025, 0.023),
        source="illustrative values; replace with measured Jet-A batch data",
    ),
    feed_line=LNGFeedLine(
        length_m=2.0,
        inner_diameter_m=12.0e-3,
        measured_heat_leak_w_per_m=10.0,
    ),
    relaxation_settings=RelaxationFlowSettings(
        relaxation_time_s=50.0e-6,
        nucleation_pressure_delay_pa=25_000.0,
    ),
    mission_fuel_masses=MissionFuelMasses(
        jet_a_kg=466.0,
        lng_kg=605.0,
        source="illustrative mission totals",
    ),
)

for point_result in result.operating_points:
    lng = point_result.lng
    jet_a = point_result.jet_a
    print(f"\n{point_result.operating_point.name}")
    if lng is not None:
        print(f"  LNG regime: {lng.regime}")
        print(
            "  LNG required diameter per orifice: "
            f"{1e3 * lng.tier2.required_orifice_diameter_m:.3f} mm"
        )
        print(f"  LNG Tier 2 exit quality: {lng.actual_exit_vapor_quality_mass:.4f}")
    if jet_a is not None:
        print(f"  Jet-A required exit diameter: {1e3 * jet_a.required_exit_diameter_m:.3f} mm")
        print(f"  Jet-A nominal cone angle: {jet_a.full_cone_angle_deg:.1f} deg")
    print(f"  warnings: {len(point_result.warnings)}")

if result.integrated_fuel_masses is not None:
    print("\nIntegrated schedule")
    print(f"  Jet-A: {result.integrated_fuel_masses.jet_a_kg:.1f} kg")
    print(f"  LNG: {result.integrated_fuel_masses.lng_kg:.1f} kg")
print(f"Study warnings: {len(result.warnings)}")