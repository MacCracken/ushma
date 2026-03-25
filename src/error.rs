//! Error types for ushma.

#[derive(Debug, thiserror::Error)]
#[non_exhaustive]
pub enum UshmaError {
    #[error("temperature must be positive: {kelvin} K (absolute zero is 0 K)")]
    InvalidTemperature { kelvin: f64 },

    #[error("pressure must be positive: {pascals} Pa")]
    InvalidPressure { pascals: f64 },

    #[error("volume must be positive: {cubic_meters} m³")]
    InvalidVolume { cubic_meters: f64 },

    #[error("negative thermal conductivity: {value} W/(m⋅K)")]
    InvalidConductivity { value: f64 },

    #[error("negative specific heat: {value} J/(kg⋅K)")]
    InvalidSpecificHeat { value: f64 },

    #[error("phase transition undefined at T={temperature} K, P={pressure} Pa")]
    UndefinedPhaseTransition { temperature: f64, pressure: f64 },

    #[error("division by zero in thermal calculation: {context}")]
    DivisionByZero { context: String },

    #[error("invalid parameter: {reason}")]
    InvalidParameter { reason: String },
}

pub type Result<T> = std::result::Result<T, UshmaError>;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_invalid_temperature() {
        let e = UshmaError::InvalidTemperature { kelvin: -10.0 };
        let msg = e.to_string();
        assert!(msg.contains("-10"));
        assert!(msg.contains("positive"));
    }

    #[test]
    fn test_invalid_pressure() {
        let e = UshmaError::InvalidPressure { pascals: -100.0 };
        assert!(e.to_string().contains("-100"));
    }

    #[test]
    fn test_invalid_volume() {
        let e = UshmaError::InvalidVolume {
            cubic_meters: -0.5,
        };
        assert!(e.to_string().contains("-0.5"));
    }

    #[test]
    fn test_invalid_conductivity() {
        let e = UshmaError::InvalidConductivity { value: -1.0 };
        assert!(e.to_string().contains("-1"));
    }

    #[test]
    fn test_invalid_specific_heat() {
        let e = UshmaError::InvalidSpecificHeat { value: -500.0 };
        assert!(e.to_string().contains("-500"));
    }

    #[test]
    fn test_phase_transition() {
        let e = UshmaError::UndefinedPhaseTransition {
            temperature: 373.15,
            pressure: 101325.0,
        };
        let msg = e.to_string();
        assert!(msg.contains("373.15"));
    }

    #[test]
    fn test_division_by_zero() {
        let e = UshmaError::DivisionByZero {
            context: "heat flux".into(),
        };
        assert!(e.to_string().contains("heat flux"));
    }

    #[test]
    fn test_invalid_parameter() {
        let e = UshmaError::InvalidParameter {
            reason: "mole count must be positive".into(),
        };
        assert!(e.to_string().contains("mole count"));
    }

    #[test]
    fn test_result_alias() {
        let ok: Result<f64> = Ok(300.0);
        assert!(ok.is_ok());
        let err: Result<f64> = Err(UshmaError::InvalidTemperature { kelvin: -1.0 });
        assert!(err.is_err());
    }

    #[test]
    fn test_error_is_debug() {
        let e = UshmaError::InvalidTemperature { kelvin: 0.0 };
        let debug = format!("{:?}", e);
        assert!(debug.contains("InvalidTemperature"));
    }
}
