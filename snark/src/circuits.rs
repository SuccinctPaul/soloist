use ark_ff::{Field, PrimeField};
use ark_relations::r1cs::{ConstraintSynthesizer, ConstraintSystemRef, SynthesisError};

#[cfg(feature = "r1cs")]
use ark_r1cs_std::prelude::*;

use std::marker::PhantomData;

#[derive(Clone)]
pub struct SimpleSha256Circuit<F: Field> {
    pub num_constraints_target: usize,
    _field: PhantomData<F>,
}

impl<F: Field> SimpleSha256Circuit<F> {
    pub fn new(num_constraints_target: usize) -> Self {
        Self {
            num_constraints_target,
            _field: PhantomData,
        }
    }
}

impl<ConstraintF: Field> ConstraintSynthesizer<ConstraintF> for SimpleSha256Circuit<ConstraintF>
where
    ConstraintF: PrimeField,
{
    #[cfg(feature = "r1cs")]
    fn generate_constraints(
        self,
        cs: ConstraintSystemRef<ConstraintF>,
    ) -> Result<(), SynthesisError> {
        for _ in 0..self.num_constraints_target {
            let _ = cs.new_witness_variable(|| Ok(ConstraintF::zero()))?;
            cs.enforce_constraint(
                ark_relations::lc!(),
                ark_relations::lc!(),
                ark_relations::lc!(),
            )?;
        }

        Ok(())
    }

    #[cfg(not(feature = "r1cs"))]
    fn generate_constraints(
        self,
        _cs: ConstraintSystemRef<ConstraintF>,
    ) -> Result<(), SynthesisError> {
        unimplemented!("This circuit requires the 'r1cs' feature to be enabled");
    }
}
