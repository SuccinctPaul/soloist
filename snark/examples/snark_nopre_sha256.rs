// usage example for 4 sub-provers in the local environment with SHA256 circuit
// RAYON_NUM_THREADS=N RUSTFLAGS='-C target-cpu=native' target-feature=+bmi2,+adx" cargo +nightly build --release --example snark_nopre_sha256 --no-default-features --features "parallel asm"
// RAYON_NUM_THREADS=32 ./snark_nopre_sha256 0/1/2/3 ../../../snark/data/4

use ark_bn254::Bn254;
use ark_ec::pairing::Pairing;
use ark_ff::UniformRand;
use ark_poly::{EvaluationDomain, GeneralEvaluationDomain};
use ark_relations::r1cs::{ConstraintSynthesizer, ConstraintSystem};
use ark_std::rand::{rngs::StdRng, SeedableRng};
use de_network::{DeMultiNet as Net, DeNet};
use merlin::Transcript;
use my_ipa::r1cs::R1CSVectors;
use my_ipa::{
    helper::{generate_r1cs_de_polynomials, generate_r1cs_pub_polynomials},
    r1cs::R1CSPubVectors,
};
use my_kzg::biv_batch_kzg::BivBatchKZG;
use my_snark::circuits::SimpleSha256Circuit;
use my_snark::snark_linear::DeSNARKLinear;
use rayon::prelude::*;
use std::path::PathBuf;
use std::time::Instant;
use structopt::StructOpt;

type MyField = <Bn254 as Pairing>::ScalarField;

#[derive(Debug, StructOpt)]
#[structopt(name = "example", about = "An example of StructOpt usage.")]
struct Opt {
    id: usize,

    #[structopt(parse(from_os_str))]
    input: PathBuf,

    nv: usize,
}

fn init() -> (usize, usize, usize) {
    let opt = Opt::from_args();
    println!("{:?}", opt);
    Net::init_from_file(opt.input.to_str().unwrap(), opt.id);
    let l = Net::n_parties();
    let sub_prover_id = Net::party_id();
    let m = 1 << opt.nv;
    (m, l, sub_prover_id)
}

fn main() {
    let (m, l, sub_prover_id) = init();
    let mut rng = StdRng::seed_from_u64(0u64);

    let challenge_r = MyField::rand(&mut rng);
    let domain_x = <GeneralEvaluationDomain<MyField> as EvaluationDomain<MyField>>::new(m).unwrap();
    let domain_y = <GeneralEvaluationDomain<MyField> as EvaluationDomain<MyField>>::new(l).unwrap();

    let x_degree = m - 1;
    let y_degree = l - 1;
    let time = Instant::now();
    let ((powers, x_srs, y_srs), v_srs) = BivBatchKZG::<Bn254>::read_or_setup(
        &mut rng,
        sub_prover_id,
        x_degree,
        y_degree,
        &domain_x,
        &domain_y,
    );
    println!("Setup time: {:?}", time.elapsed());

    let time = Instant::now();

    let circuit = SimpleSha256Circuit::<MyField>::new(m * l);

    let cs = ConstraintSystem::<MyField>::new_ref();
    circuit.clone().generate_constraints(cs.clone()).unwrap();

    assert!(cs.is_satisfied().unwrap());
    let cs_matrix = cs.to_matrices().unwrap();

    println!("Number of constraints: {:?}", cs.num_constraints());
    println!(
        "Number of variables: {:?}",
        cs.num_witness_variables() + cs.num_instance_variables()
    );

    let r1cs_vecs_all: Vec<R1CSVectors<Bn254>> = (0..l)
        .map(|sub_prover_id| {
            R1CSVectors::<Bn254>::build(sub_prover_id, m, l, challenge_r, &cs, &cs_matrix).unwrap()
        })
        .collect();
    let r1cs_de_vecs = r1cs_vecs_all[sub_prover_id].clone();
    let r1cs_de_pub_vecs: Vec<R1CSPubVectors<Bn254>> = r1cs_vecs_all
        .par_iter()
        .map(|vec| R1CSPubVectors {
            vec_x: vec.vec_x.clone(),
            vec_y: vec.vec_y.clone(),
            vec_z: vec.vec_z.clone(),
        })
        .collect();
    println!("Generate R1CS instances time: {:?}", time.elapsed());

    let time = Instant::now();
    let mut transcript: Transcript = Transcript::new(b"SHA256 R1CS inner product");
    let (sub_pub_polys, sub_wit_polys) = generate_r1cs_de_polynomials::<Bn254>(m, l, r1cs_de_vecs);
    println!("Prover {:?} starts prove", sub_prover_id);
    let proof = DeSNARKLinear::<Bn254>::de_r1cs_prove(
        sub_prover_id,
        &powers,
        &x_srs,
        &y_srs,
        &sub_wit_polys,
        &sub_pub_polys,
        &challenge_r,
        &domain_x,
        &domain_y,
        &mut transcript,
    );
    println!(
        "Prover {:?} prove total time: {:?}",
        sub_prover_id,
        time.elapsed()
    );

    if Net::am_master() {
        let proof_size = DeSNARKLinear::<Bn254>::get_proof_size(proof.as_ref().unwrap());
        println!("Proof size is {:?} bytes", proof_size);
    }

    let total_time = Instant::now();
    if Net::am_master() {
        let time = Instant::now();
        let de_pub_polys = generate_r1cs_pub_polynomials(&r1cs_de_pub_vecs);
        println!(
            "Verifier computes public polynomials time: {:?}",
            time.elapsed()
        );
        let mut transcript: Transcript = Transcript::new(b"SHA256 R1CS inner product");
        let is_valid = DeSNARKLinear::<Bn254>::r1cs_verify_no_preprocess(
            &v_srs,
            &proof.unwrap(),
            &domain_x,
            &domain_y,
            &de_pub_polys,
            &challenge_r,
            &mut transcript,
        );
        assert!(is_valid);
    }
    println!("Verify time: {:?}", total_time.elapsed());

    // Print network statistics
    let stats = Net::stats();
    println!("\n=== Network Statistics ===");
    println!("Peer ID: {}", sub_prover_id);
    println!("Bytes sent: {} bytes", stats.bytes_sent);
    println!("Bytes received: {} bytes", stats.bytes_recv);
    println!("Broadcasts: {}", stats.broadcasts);
    println!("Messages to master: {}", stats.to_master);
    println!("Messages from master: {}", stats.from_master);

    Net::deinit();
}
