package gwas

import (
	"fmt"
	"math/bits"
	"sync"
	"unsafe"

	"github.com/hhcho/sfgwas-lmm/crypto"
	"github.com/ldsec/lattigo/v2/ring"

	"math"

	"github.com/ldsec/lattigo/v2/ckks"

	"gonum.org/v1/gonum/mat"

	"go.dedis.ch/onet/v3/log"
)

type matmulPlainInnerFn func(*crypto.CryptoParams, crypto.CipherVector, crypto.PlainMatrix, int) crypto.CipherVector
type matmulInnerFn func(*crypto.CryptoParams, crypto.CipherVector, crypto.CipherMatrix, int) crypto.CipherVector

type uint128 struct {
	hi uint64
	lo uint64
}

type CipherAccV1 struct {
	acc00 []uint128
	acc01 []uint128
	acc10 []uint128
	acc11 []uint128
}

type CipherAccV2 struct {
	acc0 [][]uint128
	acc1 [][]uint128
}

type CipherVectorAccV1 []CipherAccV1

type CipherVectorAccV2 struct {
	val []CipherAccV2
	mux sync.Mutex
}

func NewCipherVectorAccV1(cryptoParams *crypto.CryptoParams, n int) CipherVectorAccV1 {
	N := cryptoParams.Params.N()
	out := make(CipherVectorAccV1, n)
	for i := range out {
		out[i].acc00 = make([]uint128, N)
		out[i].acc01 = make([]uint128, N)
		out[i].acc10 = make([]uint128, N)
		out[i].acc11 = make([]uint128, N)
	}
	return out
}

func (acc *CipherVectorAccV2) InitCipherVectorAccV2(cryptoParams *crypto.CryptoParams, n int, level int) {
	N := cryptoParams.Params.N()
	acc.val = make([]CipherAccV2, n)
	for i := range acc.val {
		acc.val[i].acc0 = make([][]uint128, level)
		acc.val[i].acc1 = make([][]uint128, level)
		for l := 0; l < level; l++ {
			acc.val[i].acc0[l] = make([]uint128, N)
			acc.val[i].acc1[l] = make([]uint128, N)
		}
	}
}

func NewCipherVectorAccV2(cryptoParams *crypto.CryptoParams, n int, level int) CipherVectorAccV2 {
	var out CipherVectorAccV2
	out.InitCipherVectorAccV2(cryptoParams, n, level)
	return out
}

func MulCoeffsAndAdd128(a, b []uint64, c []uint128) {

	var hi, lo, carry uint64

	for j := 0; j < len(a); j = j + 8 {

		x := (*[8]uint64)(unsafe.Pointer(&a[j]))
		y := (*[8]uint64)(unsafe.Pointer(&b[j]))
		z := (*[8]uint128)(unsafe.Pointer(&c[j]))

		hi, lo = bits.Mul64(x[0], y[0])
		z[0].lo, carry = bits.Add64(z[0].lo, lo, 0)
		z[0].hi += hi + carry

		hi, lo = bits.Mul64(x[1], y[1])
		z[1].lo, carry = bits.Add64(z[1].lo, lo, 0)
		z[1].hi += hi + carry

		hi, lo = bits.Mul64(x[2], y[2])
		z[2].lo, carry = bits.Add64(z[2].lo, lo, 0)
		z[2].hi += hi + carry

		hi, lo = bits.Mul64(x[3], y[3])
		z[3].lo, carry = bits.Add64(z[3].lo, lo, 0)
		z[3].hi += hi + carry

		hi, lo = bits.Mul64(x[4], y[4])
		z[4].lo, carry = bits.Add64(z[4].lo, lo, 0)
		z[4].hi += hi + carry

		hi, lo = bits.Mul64(x[5], y[5])
		z[5].lo, carry = bits.Add64(z[5].lo, lo, 0)
		z[5].hi += hi + carry

		hi, lo = bits.Mul64(x[6], y[6])
		z[6].lo, carry = bits.Add64(z[6].lo, lo, 0)
		z[6].hi += hi + carry

		hi, lo = bits.Mul64(x[7], y[7])
		z[7].lo, carry = bits.Add64(z[7].lo, lo, 0)
		z[7].hi += hi + carry
	}
}

func ReduceAndAddUint128(in []uint128, out []uint64, qInv, q uint64) {

	var hhi uint64

	for j := 0; j < len(in); j = j + 8 {

		x := (*[8]uint128)(unsafe.Pointer(&in[j]))
		y := (*[8]uint64)(unsafe.Pointer(&out[j]))

		hhi, _ = bits.Mul64(x[0].lo*qInv, q)
		y[0] += x[0].hi - hhi + q

		hhi, _ = bits.Mul64(x[1].lo*qInv, q)
		y[1] += x[1].hi - hhi + q

		hhi, _ = bits.Mul64(x[2].lo*qInv, q)
		y[2] += x[2].hi - hhi + q

		hhi, _ = bits.Mul64(x[3].lo*qInv, q)
		y[3] += x[3].hi - hhi + q

		hhi, _ = bits.Mul64(x[4].lo*qInv, q)
		y[4] += x[4].hi - hhi + q

		hhi, _ = bits.Mul64(x[5].lo*qInv, q)
		y[5] += x[5].hi - hhi + q

		hhi, _ = bits.Mul64(x[6].lo*qInv, q)
		y[6] += x[6].hi - hhi + q

		hhi, _ = bits.Mul64(x[7].lo*qInv, q)
		y[7] += x[7].hi - hhi + q
	}
}

func ModularReduceV1(cryptoParams *crypto.CryptoParams, cva CipherVectorAccV1, outScale float64) crypto.CipherVector {
	N := cryptoParams.Params.N()
	ringQ, _ := ring.NewRing(N, cryptoParams.Params.Qi())

	out := make(crypto.CipherVector, len(cva))
	for i := range out {
		ct := ckks.NewCiphertext(cryptoParams.Params, 1, 1, outScale)
		ReduceAndAddUint128(cva[i].acc00, ct.Value()[0].Coeffs[0], ringQ.MredParams[0], ringQ.Modulus[0])
		ReduceAndAddUint128(cva[i].acc01, ct.Value()[1].Coeffs[0], ringQ.MredParams[0], ringQ.Modulus[0])
		ReduceAndAddUint128(cva[i].acc10, ct.Value()[0].Coeffs[1], ringQ.MredParams[1], ringQ.Modulus[1])
		ReduceAndAddUint128(cva[i].acc11, ct.Value()[1].Coeffs[1], ringQ.MredParams[1], ringQ.Modulus[1])
		out[i] = ct
	}

	return out
}

func ModularReduceV2(cryptoParams *crypto.CryptoParams, cva CipherVectorAccV2, outScale float64) crypto.CipherVector {
	N := cryptoParams.Params.N()
	ringQ, _ := ring.NewRing(N, cryptoParams.Params.Qi())
	level := len(cva.val[0].acc0)

	out := make(crypto.CipherVector, len(cva.val))
	for i := range out {
		ct := ckks.NewCiphertext(cryptoParams.Params, 1, level-1, outScale)
		for l := 0; l < level; l++ {
			mredParams := ringQ.MredParams[l]
			qi := ringQ.Modulus[l]
			ReduceAndAddUint128(cva.val[i].acc0[l], ct.Value()[0].Coeffs[l], mredParams, qi)
			ReduceAndAddUint128(cva.val[i].acc1[l], ct.Value()[1].Coeffs[l], mredParams, qi)
		}
		err := cryptoParams.WithEvaluator(func(eval ckks.Evaluator) error {
			return eval.Reduce(ct, ct)
		})
		if err != nil {
			panic(err)
		}
		out[i] = ct
	}
	return out
}

// Multiply X and Y to add to Acc without modular reduction
func CPMultAccWithoutMRedV1(cryptoParams *crypto.CryptoParams, X crypto.CipherVector, Y crypto.PlainVector, Acc CipherVectorAccV1) {
	for i := range X {
		if X[i] != nil && Y[i] != nil {
			MulCoeffsAndAdd128(X[i].Value()[0].Coeffs[0], Y[i].Value()[0].Coeffs[0], Acc[i].acc00)
			MulCoeffsAndAdd128(X[i].Value()[1].Coeffs[0], Y[i].Value()[0].Coeffs[0], Acc[i].acc01)
			MulCoeffsAndAdd128(X[i].Value()[0].Coeffs[1], Y[i].Value()[0].Coeffs[1], Acc[i].acc10)
			MulCoeffsAndAdd128(X[i].Value()[1].Coeffs[1], Y[i].Value()[0].Coeffs[1], Acc[i].acc11)
		}
	}
}

func CPMultAccWithoutMRedV2(X crypto.CipherVector, Y crypto.PlainVector, Acc CipherVectorAccV2) {
	n := len(Acc.val)
	for i := 0; i < n; i++ {
		// Broadcasting
		xi, yi := i, i
		if len(X) == 1 {
			xi = 0
		}
		if len(Y) == 1 {
			yi = 0
		}

		if X[xi] != nil && Y[yi] != nil {
			for l := 0; l < len(Acc.val[i].acc0); l++ {
				MulCoeffsAndAdd128(X[xi].Value()[0].Coeffs[l], Y[yi].Value()[0].Coeffs[l], Acc.val[i].acc0[l])
				MulCoeffsAndAdd128(X[xi].Value()[1].Coeffs[l], Y[yi].Value()[0].Coeffs[l], Acc.val[i].acc1[l])
			}
		}
	}
}

func ToMontgomeryForm(cryptoParams *crypto.CryptoParams, pt crypto.PlainVector) {
	N := cryptoParams.Params.N()
	ringQ, _ := ring.NewRing(N, cryptoParams.Params.Qi())
	for i := range pt {
		if pt[i] != nil {
			MFormLvl(ringQ, pt[i].Level(), pt[i].Value()[0], pt[i].Value()[0])
		}
	}
}

func MFormLvl(r *ring.Ring, level int, p1, p2 *ring.Poly) {
	for i := 0; i < level+1; i++ {
		qi := r.Modulus[i]
		bredParams := r.BredParams[i]
		p1tmp, p2tmp := p1.Coeffs[i], p2.Coeffs[i]
		for j := 0; j < r.N; j = j + 8 {

			x := (*[8]uint64)(unsafe.Pointer(&p1tmp[j]))
			z := (*[8]uint64)(unsafe.Pointer(&p2tmp[j]))

			z[0] = MForm(x[0], qi, bredParams)
			z[1] = MForm(x[1], qi, bredParams)
			z[2] = MForm(x[2], qi, bredParams)
			z[3] = MForm(x[3], qi, bredParams)
			z[4] = MForm(x[4], qi, bredParams)
			z[5] = MForm(x[5], qi, bredParams)
			z[6] = MForm(x[6], qi, bredParams)
			z[7] = MForm(x[7], qi, bredParams)
		}
	}
}

func MForm(a, q uint64, u []uint64) (r uint64) {
	mhi, _ := bits.Mul64(a, u[1])
	r = -(a*u[0] + mhi) * q
	if r >= q {
		r -= q
	}
	return
}

func CMatMult1Parallel(cryptoParams *crypto.CryptoParams, A crypto.CipherMatrix, B crypto.CipherMatrix) crypto.CipherMatrix {
	m := len(B)
	s := len(A)
	slots := cryptoParams.GetSlots()
	s_ct := int((s-1)/slots) + 1
	out := crypto.CZeroMat(cryptoParams, s_ct, m)

	nproc := cryptoParams.GetNumThreads()
	fmt.Println("cmatmult1 par ", nproc)
	// Dispatcher
	jobChannels := make([]chan [2]int, nproc)
	for i := range jobChannels {
		jobChannels[i] = make(chan [2]int, 32)
	}
	go func() {
		for i := 0; i < s; i++ {
			for j := 0; j < m; j++ {
				jobChannels[(i*m+j)%nproc] <- [2]int{i, j}
			}
		}
		for _, c := range jobChannels {
			close(c)
		}
	}()

	// Workers
	var workerGroup sync.WaitGroup
	locks := make([][]sync.Mutex, m)
	for j := range locks {
		locks[j] = make([]sync.Mutex, s_ct)
	}
	for thread := 0; thread < nproc; thread++ {
		workerGroup.Add(1)
		go func(thread int, cps *crypto.CryptoParams) {
			defer workerGroup.Done()

			for pair := range jobChannels[thread] {
				i, j := pair[0], pair[1]

				ctid := i / slots
				slotid := i % slots
				prod := crypto.CMult(cps, A[i], B[j])
				ct := crypto.InnerSumAll(cps, prod)
				ct = crypto.Mask(cps, ct, slotid, false)

				locks[j][ctid].Lock()
				cps.WithEvaluator(func(eval ckks.Evaluator) error {
					eval.Add(out[j][ctid], ct, out[j][ctid])
					return nil
				})
				locks[j][ctid].Unlock()
			}
		}(thread, cryptoParams.GetThread(thread))
	}
	workerGroup.Wait()

	return out
}

func CPMatMult1Parallel(cryptoParams *crypto.CryptoParams, A crypto.PlainMatrix, B crypto.CipherMatrix) crypto.CipherMatrix {
	m := len(B)
	s := len(A)
	slots := cryptoParams.GetSlots()
	s_ct := int((s-1)/slots) + 1
	// totalVecLen := s_ct*slots
	out := crypto.CZeroMat(cryptoParams, s_ct, m)

	nproc := cryptoParams.GetNumThreads()

	// Dispatcher
	jobChannels := make([]chan [2]int, nproc)
	for i := range jobChannels {
		jobChannels[i] = make(chan [2]int, 32)
	}
	go func() {
		for i := 0; i < s; i++ {
			for j := 0; j < m; j++ {
				jobChannels[(i*m+j)%nproc] <- [2]int{i, j}
			}
		}
		for _, c := range jobChannels {
			close(c)
		}
	}()

	// Workers
	var workerGroup sync.WaitGroup
	locks := make([][]sync.Mutex, m)
	for j := range locks {
		locks[j] = make([]sync.Mutex, s_ct)
	}
	for thread := 0; thread < nproc; thread++ {
		workerGroup.Add(1)
		go func(thread int, cps *crypto.CryptoParams) {
			defer workerGroup.Done()

			for pair := range jobChannels[thread] {
				i, j := pair[0], pair[1]

				ctid := i / slots
				slotid := i % slots
				prod := crypto.CPMult(cps, B[j], A[i])
				ct := crypto.InnerSumAll(cps, prod)
				ct = crypto.Mask(cps, ct, slotid, false)

				locks[j][ctid].Lock()
				cps.WithEvaluator(func(eval ckks.Evaluator) error {
					eval.Add(out[j][ctid], ct, out[j][ctid])
					return nil
				})
				locks[j][ctid].Unlock()
			}
		}(thread, cryptoParams.GetThread(thread))
	}
	workerGroup.Wait()

	return out
}

func CMatMult5Parallel(cryptoParams *crypto.CryptoParams, A crypto.CipherMatrix, B crypto.CipherMatrix) crypto.CipherMatrix {
	n := len(A)
	m := len(B)
	s_ct := len(A[0])
	slots := cryptoParams.GetSlots()
	out := crypto.CZeroMat(cryptoParams, s_ct, m)

	nproc := cryptoParams.GetNumThreads()

	// Dispatcher
	jobChannels := make([]chan [2]int, nproc)
	for i := range jobChannels {
		jobChannels[i] = make(chan [2]int, 32)
	}
	go func() {
		for j := 0; j < n; j++ {
			for i := 0; i < m; i++ {
				jobChannels[(j*m+i)%nproc] <- [2]int{i, j}
			}
		}
		for _, c := range jobChannels {
			close(c)
		}
	}()

	// Workers
	var workerGroup sync.WaitGroup
	locks := make([]sync.Mutex, m)
	for thread := 0; thread < nproc; thread++ {
		workerGroup.Add(1)
		go func(thread int, cps *crypto.CryptoParams) {
			defer workerGroup.Done()

			for pair := range jobChannels[thread] {
				i, j := pair[0], pair[1]

				ctid := j / slots
				slotid := j % slots

				ct := crypto.Mask(cps, B[i][ctid], slotid, false)
				ct = crypto.RotateAndAdd(cps, ct, cps.GetSlots())

				res := crypto.CMult(cps, A[j], crypto.CipherVector{ct})

				locks[i].Lock()
				cps.WithEvaluator(func(eval ckks.Evaluator) error {
					for r := range res {
						eval.Add(out[i][r], res[r], out[i][r])
					}
					return nil
				})
				locks[i].Unlock()
			}
		}(thread, cryptoParams.GetThread(thread))
	}
	workerGroup.Wait()

	return out
}

func CPMatMult5Parallel(cryptoParams *crypto.CryptoParams, A crypto.PlainMatrix, B crypto.CipherMatrix) crypto.CipherMatrix {
	n := len(A)
	m := len(B)
	s_ct := len(A[0])
	slots := cryptoParams.GetSlots()
	out := crypto.CZeroMat(cryptoParams, s_ct, m)

	nproc := cryptoParams.GetNumThreads()

	// Dispatcher
	jobChannels := make([]chan [2]int, nproc)
	for i := range jobChannels {
		jobChannels[i] = make(chan [2]int, 32)
	}
	go func() {
		for j := 0; j < n; j++ {
			for i := 0; i < m; i++ {
				jobChannels[(j*m+i)%nproc] <- [2]int{i, j}
			}
		}
		for _, c := range jobChannels {
			close(c)
		}
	}()

	// Workers
	var workerGroup sync.WaitGroup
	locks := make([]sync.Mutex, m)
	for thread := 0; thread < nproc; thread++ {
		workerGroup.Add(1)
		go func(thread int, cps *crypto.CryptoParams) {
			defer workerGroup.Done()

			for pair := range jobChannels[thread] {
				i, j := pair[0], pair[1]

				ctid := j / slots
				slotid := j % slots

				ct := crypto.Mask(cps, B[i][ctid], slotid, false)
				ct = crypto.RotateAndAdd(cps, ct, cps.GetSlots())

				res := crypto.CPMult(cps, crypto.CipherVector{ct}, A[j])

				locks[i].Lock()
				cps.WithEvaluator(func(eval ckks.Evaluator) error {
					for r := range res {
						eval.Add(out[i][r], res[r], out[i][r])
					}
					return nil
				})
				locks[i].Unlock()
			}
		}(thread, cryptoParams.GetThread(thread))
	}
	workerGroup.Wait()

	return out
}
func CMatMult5(cryptoParams *crypto.CryptoParams, A crypto.CipherMatrix, B crypto.CipherMatrix) crypto.CipherMatrix {
	m := len(B)
	n := len(A) //n=C
	s_ct := len(A[0])
	slots := cryptoParams.GetSlots()
	// s_ct := int((s-1)/slots) + 1
	out := crypto.CZeroMat(cryptoParams, s_ct, m)

	for i := 0; i < m; i++ {
		for j := 0; j < n; j++ {
			ctid := int(j / slots)
			slotid := j % slots

			ct := crypto.Mask(cryptoParams, B[i][ctid], slotid, false)
			fullCt := crypto.RotateAndPlace(cryptoParams, ct, slots, 0, true)
			dupCt := make(crypto.CipherVector, s_ct)
			for ci := range dupCt {
				dupCt[ci] = fullCt
			}
			res := crypto.CMult(cryptoParams, A[j], dupCt)
			out[i] = crypto.CAdd(cryptoParams, out[i], res)
		}
	}

	return out
}

func CPMatMult5(cryptoParams *crypto.CryptoParams, A crypto.PlainMatrix, B crypto.CipherMatrix) crypto.CipherMatrix {
	n := len(A)
	m := len(B)
	s_ct := len(A[0])
	slots := cryptoParams.GetSlots()

	out := crypto.CZeroMat(cryptoParams, s_ct, m)

	for i := 0; i < m; i++ {
		for j := 0; j < n; j++ {
			ctid := j / slots
			slotid := j % slots

			ct := crypto.Mask(cryptoParams, B[i][ctid], slotid, false)
			ct = crypto.RotateAndAdd(cryptoParams, ct, slots)
			res := crypto.CPMult(cryptoParams, crypto.CipherVector{ct}, A[j])
			out[i] = crypto.CAdd(cryptoParams, out[i], res)
		}
	}

	return out
}

type BlockI8 struct {
	data [][]int8
	r    int
	c    int
}

func NewBlockI8(r, c int) BlockI8 {
	return BlockI8{
		data: make([][]int8, r),
		r:    r,
		c:    c,
	}
}

func (b BlockI8) At(i, j int) float64 {
	return float64(b.data[i][j])
}

func (b BlockI8) Dims() (int, int) {
	return b.r, b.c
}

type BlockPlainF64 struct {
	data [][]float64
	r    int
	c    int
}

func NewBlockPlainF64(r, c int) BlockPlainF64 {
	return BlockPlainF64{
		data: make([][]float64, r),
		r:    r,
		c:    c,
	}
}

func (b BlockPlainF64) At(i, j int) float64 {
	return b.data[i][j]
}

func (b BlockPlainF64) Dims() (int, int) {
	return b.r, b.c
}

type BlockF64 mat.Dense

func (b BlockF64) At(i, j int) float64 {
	return b.At(i, j)
}

type Block interface {
	At(int, int) float64
	Dims() (int, int)
}

type BlockVector []Block
type BlockMatrix []BlockVector

func ToBlockMatrix(A *mat.Dense, d int) BlockMatrix {
	r, c := A.Dims()
	br, bc := int((r-1)/d)+1, int((c-1)/d)+1

	out := make(BlockMatrix, br)
	for bi := range out {
		out[bi] = make(BlockVector, bc)
		for bj := range out[bi] {
			i1, i2 := bi*d, Min((bi+1)*d, r)
			j1, j2 := bj*d, Min((bj+1)*d, c)
			out[bi][bj] = A.Slice(i1, i2, j1, j2)
		}
	}

	return out
}

// Return if a diagonal vector exists without extracting elements
func GetDiagBool(X Block, dim int, index int) bool {
	r, c := X.Dims()
	index = Mod(index, dim) // range [0, dim-1]
	return (dim+1-r) <= index || index <= c-1
}

// index 0 is the main diagonal
// max size of Block is dim by dim and index ranges from 0 to dim-1 (mod dim)
// If given diagonal does not overlap with X (matrix might be smaller), returns false
func GetDiag(dst []float64, X Block, dim int, index int) bool {
	r, c := X.Dims()

	index = Mod(index, dim) // range [0, dim-1]

	if (dim+1-r) <= index || index <= c-1 {

		if dst == nil {
			dst = make([]float64, dim)
		} else if len(dst) < c {
			panic("destination array is not large enough")
		}

		i := Mod(-index, dim)
		for j := 0; j < len(dst); j++ {
			if i < r && j < c {
				dst[j] = X.At(i, j)
			} else {
				dst[j] = 0
			}

			i = Mod(i+1, dim)
		}

		return true
	}

	return false
}

func convertToComplex128WithRot(v []float64, nrot int) []complex128 {
	res := make([]complex128, len(v))
	for i, el := range v {
		res[Mod(i+nrot, len(res))] = complex(el, 0)
	}
	return res
}

// Return if a diagonal vector exists without extracting/encoding the vectors
func EncodeDiagBool(X BlockVector, index int, slots int) bool {
	for i := range X {
		if GetDiagBool(X[i], slots, index) {
			return true
		}
	}
	return false
}

// index specifies which diagonal to extract
// applies right-rotation by nrot positions before encoding
func EncodeDiag(cryptoParams *crypto.CryptoParams, X BlockVector, index int, nrot int, level int) (crypto.PlainVector, bool) {
	slots := cryptoParams.GetSlots()

	buf := make([]float64, slots)
	out := make(crypto.PlainVector, len(X))
	anyFlag := false

	for i := range X {
		success := GetDiag(buf, X[i], slots, index)
		if success {
			anyFlag = true
			plaintext := ckks.NewPlaintext(cryptoParams.Params, level, cryptoParams.Params.Scale())
			cryptoParams.WithEncoder(func(encoder ckks.Encoder) error {
				encoder.EncodeNTT(plaintext, convertToComplex128WithRot(buf, nrot), cryptoParams.Params.LogSlots())
				return nil
			})
			out[i] = plaintext
		} else {
			out[i] = nil
		}
	}

	return out, anyFlag
}

func EncodeDiagWithEncoder(cryptoParams *crypto.CryptoParams, X BlockVector, index int, nrot int, level int, enc ckks.Encoder) (crypto.PlainVector, bool) {
	slots := cryptoParams.GetSlots()

	buf := make([]float64, slots)
	out := make(crypto.PlainVector, len(X))
	anyFlag := false

	for i := range X {
		success := GetDiag(buf, X[i], slots, index)
		if success {
			anyFlag = true
			plaintext := ckks.NewPlaintext(cryptoParams.Params, level, cryptoParams.Params.Scale())
			enc.EncodeNTT(plaintext, convertToComplex128WithRot(buf, nrot), cryptoParams.Params.LogSlots())
			out[i] = plaintext
		} else {
			out[i] = nil
		}
	}

	return out, anyFlag
}

// Cache structure for MatMult4
// First index corresponds to row index of blocks of size slots-by-slots
// Second index of indexMap corresponds to index of diagonals (0..slots-1),
// The table maps a diag index to the encoded PlainVector
// If a given index has empty data, stored element is nil
type PlainMatrixDiagCache [][]crypto.PlainVector

func MatSATASClear(gfs *GenoFileStream, invStdDev []float64) *mat.Dense {
	gfs.Reset()
	numRows := int(gfs.NumRows())
	numCols := int(gfs.NumCols())
	A := mat.NewDense(numRows, numCols, nil)
	counter := 0
	for i := 0; i < numRows; i++ {
		row := gfs.NextRow()
		rowFloat := row
		for j := 0; j < len(row); j++ {
			rowFloat[j] = rowFloat[j] * invStdDev[j]
		}
		A.SetRow(counter, rowFloat)
		counter++
		if counter%10000 == 0 {
			fmt.Println("ATA", counter)
		}
	}
	ATA := mat.NewDense(numCols, numCols, nil)
	ATA.Mul(A.T(), A)
	return ATA
}
func MatTSClearSimple(gfs *GenoFileStream, B *mat.Dense, invStdDev []float64) *mat.Dense {
	gfs.Reset()
	numRows := int(gfs.NumRows())
	numCols := int(gfs.NumCols())
	A := mat.NewDense(numRows, numCols, nil)
	counter := 0
	for i := 0; i < numRows; i++ {
		row := gfs.NextRow()
		rowFloat := row
		for j := 0; j < len(row); j++ {
			rowFloat[j] = rowFloat[j] * invStdDev[j]
		}
		A.SetRow(counter, rowFloat)
		counter++
	}
	_, c := B.Dims()
	ATB := mat.NewDense(numCols, c, nil)
	ATB.Mul(A.T(), B)
	return ATB
}

func MatTClear(gfs *GenoFileStream, B *mat.Dense) *mat.Dense {
	gfs.Reset()
	aa, c := B.Dims()
	numRows := int(gfs.NumRowsToKeep())
	numCols := int(gfs.NumCols())
	blockSize := 20 // need to set
	log.LLvl1("numRows", numRows, "numCols", numCols, "B dims", aa, c)
	if gfs.numRows < 20 {
		panic("num rows less than blocksize")
	}
	blocks := ((numRows - 1) / blockSize) + 1
	out := mat.NewDense(numCols, c, nil)
	for b := 0; b < blocks; b++ {
		realSize := blockSize
		if b == blocks-1 && numRows%blockSize != 0 {
			realSize = numRows % blockSize
		}

		block := mat.NewDense(realSize, numCols, nil)
		for i := 0; i < realSize; i++ {
			row := gfs.NextRow()
			block.SetRow(i, row)
		}
		var slicedB mat.Matrix
		if b == blocks-1 {
			slicedB = B.Slice(b*blockSize, realSize+b*blockSize, 0, c)
		} else {
			slicedB = B.Slice(b*blockSize, (b+1)*blockSize, 0, c)
		}
		newRes := mat.NewDense(numCols, c, nil)
		newRes.Mul(block.T(), slicedB)
		out.Add(out, newRes)
	}

	return out
}

func MatMult4StreamNProc(cryptoParams *crypto.CryptoParams, A crypto.CipherMatrix, gfs *GenoFileStream, maxLevel int, computeSquaredSum bool, nproc int) (crypto.CipherMatrix, []float64) {
	gfs.Reset() // Reset to beginning of file just in case

	nrow, ncol := gfs.NumRowsToKeep(), gfs.NumColsToKeep()

	s := len(A)
	outScale := A[0][0].Scale() * cryptoParams.Params.Scale()
	slots := cryptoParams.GetSlots()
	d := int(math.Ceil(math.Sqrt(float64(slots))))
	m_ct := ((ncol - 1) / uint64(slots)) + 1
	numBlockRows := ((nrow - 1) / uint64(slots)) + 1

	A, _ = crypto.FlattenLevels(cryptoParams, A)
	if A[0][0].Level() > maxLevel {
		A = crypto.DropLevel(cryptoParams, A, maxLevel)
	}

	accCache := make([][]CipherVectorAccV2, s)
	accCacheMux := make([][]sync.Mutex, s)
	for i := range accCache {
		accCache[i] = make([]CipherVectorAccV2, d) // Cache each of the sqrt(slots) groups, initialize later on-the-fly
		accCacheMux[i] = make([]sync.Mutex, d)
	}

	rotCache := make(crypto.CipherMatrix, s)
	for i := range rotCache {
		rotCache[i] = make(crypto.CipherVector, d)
	}

	var sqSum []float64
	if computeSquaredSum {
		sqSum = make([]float64, ncol)
	}

	for bi := 0; bi < int(numBlockRows); bi++ {
		BSlice := make([]BlockPlainF64, m_ct)
		nr := Min((bi+1)*slots, int(nrow)) - bi*slots
		for ri := 0; ri < nr; ri++ {

			// Read one row from file
			row := gfs.NextRow()

			for rj := range row {
				if computeSquaredSum {
					sqSum[rj] += float64(row[rj] * row[rj])
				}
			}

			// Add slice to each block matrix
			for bj := range BSlice {
				j1 := bj * slots
				j2 := Min((bj+1)*slots, int(ncol))
				nc := j2 - j1
				if ri == 0 {
					BSlice[bj] = NewBlockPlainF64(nr, nc)
				}
				BSlice[bj].data[ri] = row[j1:j2]
			}
		}

		blockVec := make(BlockVector, m_ct)
		for bj := range blockVec {
			blockVec[bj] = Block(BSlice[bj])
		}

		// Pre-collect active baby/giant indices
		babyTable := make([]bool, d)
		giantTable := make([]bool, d)
		shiftTable := make([]bool, slots)
		for shift := 0; shift < slots; shift++ {
			if EncodeDiagBool(blockVec, -shift, slots) {
				baby, giant := shift%d, shift/d
				babyTable[baby] = true
				giantTable[giant] = true
				shiftTable[shift] = true
			}
		}

		// Dispatcher
		jobChannels := make([]chan int, nproc)
		for i := range jobChannels {
			jobChannels[i] = make(chan int, 64)
		}
		go func() {
			index := 0
			for baby, flag := range babyTable {
				if flag {
					jobChannels[index%nproc] <- baby
					index++
				}
			}
			for _, c := range jobChannels {
				close(c)
			}
		}()

		// Workers
		var workerGroup sync.WaitGroup
		Aslice := make(crypto.CipherVector, len(A))
		for i := range A {
			Aslice[i] = A[i][bi]
		}
		for thread := 0; thread < nproc; thread++ {
			workerGroup.Add(1)
			go func(thread int) {
				defer workerGroup.Done()

				eva := ckks.NewEvaluator(cryptoParams.Params, ckks.EvaluationKey{Rlk: cryptoParams.Rlk, Rtks: cryptoParams.RotKs})

				for baby := range jobChannels[thread] {
					for i := range A {
						rotCache[i][baby] = crypto.RotateRightWithEvaluator(cryptoParams, Aslice[i], -baby, eva)
					}
				}
			}(thread)
		}
		workerGroup.Wait()

		for giant, flag := range giantTable {
			if flag {
				for i := range A {
					if accCache[i][giant].val == nil {
						accCache[i][giant] = NewCipherVectorAccV2(cryptoParams, int(m_ct), maxLevel)
					}
				}
			}
		}

		// Extract and encode diagonal vectors
		shiftChannels := make([]chan int, nproc)
		for i := range shiftChannels {
			shiftChannels[i] = make(chan int, 128)
		}

		go func() {
			index := 0
			for shift, flag := range shiftTable {
				if flag {
					shiftChannels[index%nproc] <- shift
					index++
				}
			}
			for _, c := range shiftChannels {
				close(c)
			}
		}()

		for thread := 0; thread < nproc; thread++ {
			workerGroup.Add(1)
			go func(thread int) {
				defer workerGroup.Done()

				enc := ckks.NewEncoderBig(cryptoParams.Params, cryptoParams.GetPrec())

				for shift := range shiftChannels[thread] {
					baby, giant := shift%d, shift/d

					plainVec, _ := EncodeDiagWithEncoder(cryptoParams, blockVec, -shift, d*giant, maxLevel, enc)

					ToMontgomeryForm(cryptoParams, plainVec)

					for i := range A {
						accCacheMux[i][giant].Lock()
						CPMultAccWithoutMRedV2(crypto.CipherVector{rotCache[i][baby]}, plainVec, accCache[i][giant])
						accCacheMux[i][giant].Unlock()
					}
				}
			}(thread)
		}
		workerGroup.Wait()
	}

	out := crypto.CZeroMat(cryptoParams, int(m_ct), s)
	for i := range out {
		jobChannels := make([]chan int, nproc)
		for j := range jobChannels {
			jobChannels[j] = make(chan int, 32)
		}

		go func() {
			for l := range accCache[i] {
				if accCache[i][l].val != nil {
					jobChannels[l%nproc] <- l
				}
			}
			for _, c := range jobChannels {
				close(c)
			}
		}()

		aggChannel := make(chan crypto.CipherVector, 8)

		var wg sync.WaitGroup
		for thread := 0; thread < nproc; thread++ {
			wg.Add(1)
			go func(thread int) {
				defer wg.Done()

				eva := ckks.NewEvaluator(cryptoParams.Params, ckks.EvaluationKey{Rlk: cryptoParams.Rlk, Rtks: cryptoParams.RotKs})

				for l := range jobChannels[thread] {
					cv := ModularReduceV2(cryptoParams, accCache[i][l], outScale)

					if l > 0 { // Giant step alignment
						for j := range cv {
							cv[j] = crypto.RotateRightWithEvaluator(cryptoParams, cv[j], -l*d, eva)
						}
					}

					aggChannel <- cv
				}
			}(thread)
		}

		var aggGroup sync.WaitGroup
		aggGroup.Add(1)
		go func() {
			defer aggGroup.Done()

			eva := ckks.NewEvaluator(cryptoParams.Params, ckks.EvaluationKey{Rlk: cryptoParams.Rlk, Rtks: cryptoParams.RotKs})

			for cv := range aggChannel {
				for j := range cv {
					eva.Add(out[i][j], cv[j], out[i][j])
				}
			}
		}()

		wg.Wait()
		close(aggChannel)
		aggGroup.Wait()
	}

	return out, sqSum
}
func MatMult4TransformB(cryptoParams *crypto.CryptoParams, B *mat.Dense) PlainMatrixDiagCache {
	slots := cryptoParams.GetSlots()
	d := int(math.Ceil(math.Sqrt(float64(slots))))
	blockB := ToBlockMatrix(B, slots)

	cache := make(PlainMatrixDiagCache, len(blockB))

	for bi := range blockB {

		cache[bi] = make([]crypto.PlainVector, slots)

		for shift := 0; shift < slots; shift++ {
			giant := int(shift / d)
			plainVec, flag := EncodeDiag(cryptoParams, blockB[bi], -shift, d*giant, cryptoParams.Params.MaxLevel())
			if !flag {
				cache[bi][shift] = nil
			} else {
				ToMontgomeryForm(cryptoParams, plainVec)
				cache[bi][shift] = plainVec
			}
		}
	}

	return cache
}

// Generalized to levels >= 2
func CPMatMult4V2CachedB(cryptoParams *crypto.CryptoParams, A crypto.CipherMatrix, maxLevel int, CachedB PlainMatrixDiagCache) crypto.CipherMatrix {
	s := len(A)
	slots := cryptoParams.GetSlots()
	d := int(math.Ceil(math.Sqrt(float64(slots))))

	if A[0][0].Level() > maxLevel {
		A = crypto.DropLevel(cryptoParams, A, maxLevel)
	}

	out := make(crypto.CipherMatrix, s)
	outScale := A[0][0].Scale() * cryptoParams.Params.Scale()

	for i := range A {

		accCache := make([]CipherVectorAccV2, d) // Cache each of the sqrt(slots) groups

		for bi := range CachedB {

			rotCache := make(crypto.CipherVector, d)

			for shift := 0; shift < slots; shift++ {
				if CachedB[bi][shift] == nil {
					continue
				}

				baby, giant := shift%d, int(shift/d)
				plainVec := CachedB[bi][shift]

				if rotCache[baby] == nil {
					rotCache[baby] = crypto.RotateRight(cryptoParams, A[i][bi], -baby)
				}

				cipherVec := make(crypto.CipherVector, len(plainVec))
				for j := range cipherVec {
					cipherVec[j] = rotCache[baby]
				}

				if accCache[giant].val == nil {
					accCache[giant] = NewCipherVectorAccV2(cryptoParams, len(plainVec), maxLevel)
				}

				CPMultAccWithoutMRedV2(cipherVec, plainVec, accCache[giant])
			}
		}

		for l := range accCache {
			if accCache[l].val != nil {

				cv := ModularReduceV2(cryptoParams, accCache[l], outScale)
				if l > 0 { // Giant step alignment
					for j := range cv {
						cv[j] = crypto.RotateRight(cryptoParams, cv[j], -l*d)
					}
				}

				if out[i] == nil {
					out[i] = cv
				} else {
					out[i] = crypto.CAdd(cryptoParams, out[i], cv)
				}
			}
		}
	}
	return out
}
