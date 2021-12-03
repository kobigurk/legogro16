use ark_ec::{AffineCurve, PairingEngine, ProjectiveCurve};
use ark_ff::Zero;
use ark_std::marker::PhantomData;
use ark_std::ops::{AddAssign, Mul};
use ark_std::vec;
use ark_std::vec::Vec;

/// CoeffPos: A struct to help build sparse matrices.
#[derive(Clone, Debug)]
pub struct CoeffPos<T> {
    val: T,
    pos: usize,
}

// a column is a vector of CoeffPos-s
type Col<T> = Vec<CoeffPos<T>>;

/* TODO: One could consider a cache-friendlier implementation for the 2-row case*/

/// Column-Major Sparse Matrix
#[derive(Clone, Debug)]
pub struct SparseMatrix<T> {
    cols: Vec<Col<T>>, // a vector of columns
    pub nr: usize,
    pub nc: usize,
}

impl<T: Copy> SparseMatrix<T> {
    // NB: Given column by column
    pub fn new(nr: usize, nc: usize) -> SparseMatrix<T> {
        SparseMatrix {
            cols: vec![vec![]; nc],
            nr,
            nc,
        }
    }

    pub fn insert_val(&mut self, r: usize, c: usize, v: &T) {
        let coeff_pos = CoeffPos { pos: r, val: *v };
        self.cols[c].push(coeff_pos);
    }

    // insert a continuous sequence of values at row r starting from c_offset
    pub fn insert_row_slice(&mut self, r: usize, c_offset: usize, vs: &[T]) {
        // NB: could be improved in efficiency by first extending the vector
        for (i, x) in vs.iter().enumerate() {
            self.insert_val(r, c_offset + i, x);
        }
    }

    pub fn get_col(&self, c: usize) -> &Col<T> {
        &self.cols[c]
    }
}

pub struct SparseLinAlgebra<PE: PairingEngine> {
    pairing_engine_type: PhantomData<PE>,
}

impl<PE: PairingEngine> SparseLinAlgebra<PE> {
    // this is basically a multi-exp
    pub fn sparse_inner_product(v: &Vec<PE::Fr>, w: &Col<PE::G1Affine>) -> PE::G1Affine {
        let mut res: PE::G1Projective = PE::G1Projective::zero();
        for coeffpos in w {
            let g = coeffpos.val;
            let i = coeffpos.pos;
            // XXX: Should this be optimized for special cases
            //         (e.g. 0 or 1) or is this already in .mul?
            let tmp = g.mul(v[i]);

            res.add_assign(&tmp);
        }
        res.into_affine()
    }

    pub fn sparse_vector_matrix_mult(
        v: &Vec<PE::Fr>,
        m: &SparseMatrix<PE::G1Affine>,
        t: usize,
    ) -> Vec<PE::G1Affine> {
        // the result should contain every column of m multiplied by v
        let mut res: Vec<PE::G1Affine> = Vec::with_capacity(t);
        for c in 0..m.nc {
            res.push(Self::sparse_inner_product(&v, &m.get_col(c)));
        }
        res
    }
}

pub fn inner_product<PE: PairingEngine>(v: &[PE::Fr], w: &[PE::G1Affine]) -> PE::G1Affine {
    assert_eq!(v.len(), w.len());
    let mut res: PE::G1Projective = PE::G1Projective::zero();
    for i in 0..v.len() {
        let tmp = w[i].mul(v[i]);
        res.add_assign(&tmp);
    }
    res.into_affine()
}

pub fn scalar_vector_mult<PE: PairingEngine>(a: &PE::Fr, v: &[PE::Fr], l: usize) -> Vec<PE::Fr> {
    let mut res: Vec<PE::Fr> = Vec::with_capacity(l);
    for i in 0..v.len() {
        let x: PE::Fr = a.mul(&v[i]);
        res.push(x);
    }
    res
}
