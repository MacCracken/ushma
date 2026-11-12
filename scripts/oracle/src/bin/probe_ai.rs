//! `probe_ai` — the request the Rust oracle's `ai` client sends, for M0 probe 7 (port plan §4 step 9; compared
//! with what sandhi sends, `scripts/probes/p7_sandhi_post.cyr`). It calls ushma 1.4.0's
//! `DaimonClient::register_agent` once with an API key against the loopback stub
//! (`scripts/probes/p7_daimon_stub.py`, which records the raw request bytes) and prints the outcome. 1.4.0 posts to
//! `/v1/agents/register` and reads a string `agent_id` (the P8 oracle defect), so the stub's 404 there is the
//! expected outcome: what is compared is the request.
//!
//! ```sh
//! python3 -I -B scripts/probes/p7_daimon_stub.py 18731 build/p7-rust.log &
//! cargo run --release --locked --manifest-path scripts/oracle/Cargo.toml --bin probe_ai -- http://127.0.0.1:18731
//! ```

use std::process::ExitCode;

use ushma::ai::{DaimonClient, DaimonConfig};

fn main() -> ExitCode {
    let endpoint = std::env::args()
        .nth(1)
        .unwrap_or_else(|| "http://127.0.0.1:18731".to_string());
    let config = DaimonConfig {
        endpoint,
        api_key: Some("probe-key".to_string()),
    };
    let client = match DaimonClient::new(config) {
        Ok(c) => c,
        Err(e) => {
            eprintln!("probe_ai: {e}");
            return ExitCode::FAILURE;
        }
    };
    let runtime = match tokio::runtime::Builder::new_current_thread()
        .enable_all()
        .build()
    {
        Ok(r) => r,
        Err(e) => {
            eprintln!("probe_ai: runtime: {e}");
            return ExitCode::FAILURE;
        }
    };
    match runtime.block_on(client.register_agent()) {
        Ok(id) => println!("register_agent: Ok({id:?})"),
        Err(e) => println!("register_agent: Err({e})"),
    }
    ExitCode::SUCCESS
}
