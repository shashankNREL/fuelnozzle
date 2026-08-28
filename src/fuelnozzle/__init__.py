"""Reduced-order dual-fuel nozzle design tools."""

from fuelnozzle.feed import LNGFeedLine, LNGFeedLineResult, solve_lng_feed_line
from fuelnozzle.jet_a import (
	JetAProperties,
	JetAPropertyTable,
	PressureSwirlGeometry,
	PressureSwirlResult,
	solve_jet_a_pressure_swirl,
)
from fuelnozzle.lng import (
	EquilibriumFlowSettings,
	FlashLocation,
	LNGNozzleGeometry,
	RelaxationFlowSettings,
	Tier0FlashScreen,
	Tier1EquilibriumFlow,
	Tier2RelaxationFlow,
	screen_lng_flash,
	solve_lng_equilibrium_flow,
	solve_lng_relaxation_flow,
)
from fuelnozzle.models import LNGComposition, ModelWarning, ThermodynamicState
from fuelnozzle.operating import OperatingPoint, PressureBudget, fuel_pressure_budget
from fuelnozzle.properties import CoolPropLNGProvider, PropertyCalculationError
from fuelnozzle.spray import (
	CFDSprayBoundary,
	FlashSprayCalibration,
	FlashSprayRegime,
	FlashSpraySettings,
	Tier3FlashSpray,
	solve_lng_flash_spray,
)
from fuelnozzle.study import (
	IntegratedFuelMasses,
	MissionFuelMasses,
	NozzleStudyResult,
	OperatingPointStudyResult,
	StudySettings,
	run_nozzle_study,
)

__version__ = "0.1.0"

__all__ = [
	"CoolPropLNGProvider",
	"CFDSprayBoundary",
	"EquilibriumFlowSettings",
	"FlashLocation",
	"FlashSprayCalibration",
	"FlashSprayRegime",
	"FlashSpraySettings",
	"IntegratedFuelMasses",
	"JetAProperties",
	"JetAPropertyTable",
	"LNGFeedLine",
	"LNGFeedLineResult",
	"LNGComposition",
	"LNGNozzleGeometry",
	"ModelWarning",
	"MissionFuelMasses",
	"NozzleStudyResult",
	"OperatingPoint",
	"OperatingPointStudyResult",
	"PressureBudget",
	"PressureSwirlGeometry",
	"PressureSwirlResult",
	"PropertyCalculationError",
	"RelaxationFlowSettings",
	"Tier0FlashScreen",
	"Tier1EquilibriumFlow",
	"Tier2RelaxationFlow",
	"Tier3FlashSpray",
	"ThermodynamicState",
	"StudySettings",
	"fuel_pressure_budget",
	"screen_lng_flash",
	"run_nozzle_study",
	"solve_lng_feed_line",
	"solve_lng_equilibrium_flow",
	"solve_lng_flash_spray",
	"solve_lng_relaxation_flow",
	"solve_jet_a_pressure_swirl",
]
