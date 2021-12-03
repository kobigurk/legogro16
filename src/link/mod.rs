#[macro_use]

mod matrix;
mod snark;

pub use matrix::*;
pub use snark::*;

#[cfg(test)]
mod test {
    use super::{PESubspaceSnark, SparseMatrix, SubspaceSnark, PP};
    use ark_bls12_381::{Bls12_381, Fr, G1Affine, G1Projective, G2Affine, G2Projective};
    use ark_ec::{AffineCurve, ProjectiveCurve};
    use ark_ff::{One, PrimeField, UniformRand, Zero};
    use ark_std::rand::{rngs::StdRng, SeedableRng};

    #[test]
    fn test_basic() {
        let mut rng = StdRng::seed_from_u64(0u64);
        let g1 = G1Projective::rand(&mut rng).into_affine();
        let g2 = G2Projective::rand(&mut rng).into_affine();

        let mut pp = PP::<G1Affine, G2Affine> { l: 1, t: 2, g1, g2 };

        let mut m = SparseMatrix::new(1, 2);
        m.insert_row_slice(0, 0, &vec![g1, g1]);

        let x: Vec<Fr> = vec![Fr::one(), Fr::zero()];

        let x_bad: Vec<Fr> = vec![Fr::one(), Fr::one()];

        let y: Vec<G1Affine> = vec![g1];

        let (ek, vk) = PESubspaceSnark::<Bls12_381>::keygen(&mut rng, &pp, m);

        let pi = PESubspaceSnark::<Bls12_381>::prove(&mut pp, &ek, &x);
        let pi_bad = PESubspaceSnark::<Bls12_381>::prove(&mut pp, &ek, &x_bad);

        let b = PESubspaceSnark::<Bls12_381>::verify(&pp, &vk, &y, &pi);
        let b_bad = PESubspaceSnark::<Bls12_381>::verify(&pp, &vk, &y, &pi_bad);
        assert!(b);
        assert!(!b_bad);
    }

    #[test]
    fn test_same_value_different_bases() {
        let mut rng = StdRng::seed_from_u64(0u64);
        let g1 = G1Projective::rand(&mut rng).into_affine();
        let g2 = G2Projective::rand(&mut rng).into_affine();

        let mut pp = PP::<G1Affine, G2Affine> { l: 2, t: 3, g1, g2 };

        let bases1 = [G1Projective::rand(&mut rng), G1Projective::rand(&mut rng)]
            .iter()
            .map(|p| p.into_affine())
            .collect::<Vec<_>>();
        let bases2 = [G1Projective::rand(&mut rng), G1Projective::rand(&mut rng)]
            .iter()
            .map(|p| p.into_affine())
            .collect::<Vec<_>>();
        let mut m = SparseMatrix::new(2, 3);
        m.insert_row_slice(0, 0, &vec![bases1[0]]);
        m.insert_row_slice(0, 2, &vec![bases1[1]]);
        m.insert_row_slice(1, 1, &vec![bases2[0], bases2[1]]);

        let x: Vec<Fr> = vec![Fr::rand(&mut rng), Fr::rand(&mut rng), Fr::rand(&mut rng)];

        let y: Vec<G1Affine> = vec![
            bases1[0].into_projective().mul(x[0].into_repr()) + bases1[1].mul(x[2].into_repr()),
            bases2[0].into_projective().mul(x[1].into_repr()) + bases2[1].mul(x[2].into_repr()),
        ]
        .into_iter()
        .map(|p| p.into_affine())
        .collect::<Vec<_>>();

        let (ek, vk) = PESubspaceSnark::<Bls12_381>::keygen(&mut rng, &pp, m);

        let pi = PESubspaceSnark::<Bls12_381>::prove(&mut pp, &ek, &x);

        let b = PESubspaceSnark::<Bls12_381>::verify(&pp, &vk, &y, &pi);
        assert!(b);
    }
}
