use std::process::{ExitCode, Termination};

#[derive(thiserror::Error, Debug, Clone)]
pub enum AppError {
    // #[error("Internal error.")]
    // Internal,
}

impl Termination for AppError {
    fn report(self) -> ExitCode {
        match self {
            // Internal => ExitCode::from(1),
            //   Other => ExitCode::from(255),
        }
    }
}
