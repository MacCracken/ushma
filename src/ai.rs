//! AI integration — daimon/hoosh client for ushma.

use serde::{Deserialize, Serialize};

use crate::error::{Result, UshmaError};

#[derive(Clone, Serialize, Deserialize)]
pub struct DaimonConfig {
    pub endpoint: String,
    pub api_key: Option<String>,
}

/// Redacts `api_key` so it cannot leak through `{:?}` or tracing fields.
impl std::fmt::Debug for DaimonConfig {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("DaimonConfig")
            .field("endpoint", &self.endpoint)
            .field("api_key", &self.api_key.as_ref().map(|_| "<redacted>"))
            .finish()
    }
}

impl Default for DaimonConfig {
    fn default() -> Self {
        Self {
            endpoint: "http://localhost:8090".into(),
            api_key: None,
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct HooshConfig {
    pub endpoint: String,
}

impl Default for HooshConfig {
    fn default() -> Self {
        Self {
            endpoint: "http://localhost:8088".into(),
        }
    }
}

#[derive(Debug)]
pub struct DaimonClient {
    config: DaimonConfig,
    client: reqwest::Client,
}

impl DaimonClient {
    pub fn new(config: DaimonConfig) -> Result<Self> {
        let client = reqwest::Client::builder()
            .timeout(std::time::Duration::from_secs(30))
            .build()
            .map_err(|e| UshmaError::InvalidParameter {
                reason: format!("failed to build HTTP client: {e}"),
            })?;
        Ok(Self { config, client })
    }

    /// Register ushma as an agent with daimon; returns the assigned agent id.
    #[tracing::instrument(level = "debug", skip(self), fields(endpoint = %self.config.endpoint))]
    pub async fn register_agent(&self) -> Result<String> {
        let body = serde_json::json!({
            "name": "ushma",
            "capabilities": ["heat_transfer", "thermodynamics", "entropy", "thermal_materials"],
        });
        let mut req = self
            .client
            .post(format!(
                "{}/v1/agents/register",
                self.config.endpoint.trim_end_matches('/')
            ))
            .json(&body);
        if let Some(ref key) = self.config.api_key {
            req = req.bearer_auth(key);
        }
        let resp = req
            .send()
            .await
            .and_then(reqwest::Response::error_for_status)
            .map_err(|e| UshmaError::InvalidParameter {
                reason: format!("registration request failed: {e}"),
            })?;
        let data: serde_json::Value =
            resp.json()
                .await
                .map_err(|e| UshmaError::InvalidParameter {
                    reason: format!("invalid registration response: {e}"),
                })?;
        let agent_id = data["agent_id"]
            .as_str()
            .ok_or_else(|| UshmaError::InvalidParameter {
                reason: "registration response missing agent_id field".into(),
            })?;
        Ok(agent_id.to_string())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_default_config() {
        let c = DaimonConfig::default();
        assert_eq!(c.endpoint, "http://localhost:8090");
        assert!(c.api_key.is_none());
    }

    #[test]
    fn test_hoosh_default() {
        let c = HooshConfig::default();
        assert_eq!(c.endpoint, "http://localhost:8088");
    }

    #[test]
    fn test_daimon_client_new() {
        let config = DaimonConfig::default();
        let client = DaimonClient::new(config);
        assert!(client.is_ok());
    }

    #[test]
    fn test_daimon_config_with_api_key() {
        let c = DaimonConfig {
            endpoint: "https://custom.host:9090".into(),
            api_key: Some("test-key".into()),
        };
        assert_eq!(c.api_key.as_deref(), Some("test-key"));
    }

    #[test]
    fn test_daimon_config_serde_roundtrip() {
        let c = DaimonConfig {
            endpoint: "http://example.com".into(),
            api_key: Some("key123".into()),
        };
        let json = serde_json::to_string(&c).unwrap();
        let back: DaimonConfig = serde_json::from_str(&json).unwrap();
        assert_eq!(back.endpoint, c.endpoint);
        assert_eq!(back.api_key, c.api_key);
    }

    #[test]
    fn test_hoosh_config_serde_roundtrip() {
        let c = HooshConfig::default();
        let json = serde_json::to_string(&c).unwrap();
        let back: HooshConfig = serde_json::from_str(&json).unwrap();
        assert_eq!(back.endpoint, c.endpoint);
    }

    #[test]
    fn test_daimon_config_debug_redacts_api_key() {
        let c = DaimonConfig {
            endpoint: "http://example.com".into(),
            api_key: Some("super-secret".into()),
        };
        let dbg = format!("{c:?}");
        assert!(!dbg.contains("super-secret"));
        assert!(dbg.contains("<redacted>"));
        let client = DaimonClient::new(c).unwrap();
        assert!(!format!("{client:?}").contains("super-secret"));
    }
}
