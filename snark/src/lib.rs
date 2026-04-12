#![deny(warnings, future_incompatible, nonstandard_style)]
#![allow(unused_imports)]

pub mod circuits;
pub mod gadgets_and_tests;
pub mod indexer;
pub mod prover_nopre;
pub mod prover_pre;
pub mod snark_linear;
pub mod snark_log;

mod serialize;

pub(crate) use serialize::{impl_serde_for_ark_serde_checked, impl_serde_for_ark_serde_unchecked};
