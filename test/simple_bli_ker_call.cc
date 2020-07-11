#include <iostream>
#include <blis/blis.h>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>

extern "C"
void execute_kernel(dgemm_ukr_ft krn,
                    unsigned k,
                    double *c1,
                    double *d1,
                    double *d2,
                    double *c2,
                    double *d3,
                    unsigned ld0, unsigned ld1,
                    auxinfo_t *aux, cntx_t *cntx)
{ krn(k, c1, d1, d2, c2, d3, ld0, ld1, aux, cntx); }

int main(const int argc, const char *argv[])
{
    // these two doesn't conflict.
    using namespace std;
    using namespace boost::numeric::ublas;

    // blis preparations.
    auto *cntx = bli_gks_query_cntx();
    auto *ugemm = (dgemm_ukr_ft)bli_cntx_get_l3_nat_ukr_dt(BLIS_DOUBLE, BLIS_GEMM_UKR, cntx);
    // queries block size.
    inc_t mr = bli_cntx_get_blksz_def_dt(BLIS_DOUBLE, BLIS_MR, cntx);
    inc_t nr = bli_cntx_get_blksz_def_dt(BLIS_DOUBLE, BLIS_NR, cntx);
    cout << mr << '\t' << nr << '\n';
    // writes address of queried microkernel and the most common haswell 6x8 kernel.
    // cout << (void *)ugemm << '\t' << (void *)&bli_dgemm_haswell_asm_6x8 << endl;

    // prepares input and output container.
    matrix<double, column_major> matA(mr, 8);
    matrix<double, column_major> matB(nr, 8);
    matrix<double, column_major> matC(mr, nr);
    for (unsigned j = 0; j < 8; ++j) {
        for (unsigned i = 0; i < mr; ++i)
            matA(i, j) = i * 8 + j;
        for (unsigned i = 0; i < nr; ++i)
            matB(i, j) = i + j * 8;
    }
    for (unsigned j = 0; j < nr; ++j)
        for (unsigned i = 0; i < mr; ++i)
            matC(i, j) = 0.0;
    double *ptrA = &(matA(0, 0));
    double *ptrB = &(matB(0, 0));
    double *ptrC = &(matC(0, 0));
    double alpha = 1.0, beta = 0.0;
    cout << matA << endl;
    cout << matB << endl;
    for (unsigned i = 0; i < mr*8; ++i)
        cout << ptrA[i] << ' ';
    cout << '\n';
    for (unsigned i = 0; i < nr*8; ++i)
        cout << ptrB[i] << ' ';
    cout << '\n';

    // auxilliary information.
    auxinfo_t def_data;
    bli_auxinfo_set_is_a(1, &def_data);
    bli_auxinfo_set_is_b(1, &def_data);
    // bli_auxinfo_set_ps_a(mr*8, &def_data); // 2nd idx stride, maybe...
    // bli_auxinfo_set_ps_b(nr*8, &def_data); // 2nd idx stride, maybe...
    // bli_auxinfo_set_schema_a(BLIS_PACKED_COL_PANELS, &def_data);
    // bli_auxinfo_set_schema_b(BLIS_PACKED_ROW_PANELS, &def_data);
    bli_auxinfo_set_next_ab(ptrA, ptrB, &def_data);

    // Call kernel (in form of c function pointer) from c++:
    // ugemm(8, &alpha, ptrA, ptrB, &beta, ptrC, column_skip=1, mr, def_data, cntx);
    execute_kernel(ugemm, 8, &alpha, ptrA, ptrB, &beta, ptrC, 1, mr, &def_data, cntx);
    // usually it's calling:
    // bli_dgemm_haswell_asm_6x8(8, &alpha, ptrA, ptrB, &beta, ptrC, 1, mr, &def_data, cntx);
    // equivalent to calling:
    // bli_dgemm(BLIS_NO_TRANSPOSE, BLIS_TRANSPOSE, mr, nr, 8,
    //           &alpha, ptrA, 1, mr, ptrB, 1, nr, &beta, ptrC, 1, mr);

    // row-col formatted output of C.
    for (unsigned i = 0; i < mr; i++) {
        for (unsigned j = 0; j < nr; j++)
            cout << matC(i, j) << '\t';
        cout << endl;
    }
    return 0;
}