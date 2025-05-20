package lmm

import (
	"errors"
	"fmt"
	"math"
	"os"
	"runtime"
	"sort"
	"sync"

	"strconv"
	"time"

	"go.dedis.ch/onet/v3/log"

	mpc_core "github.com/hhcho/mpc-core"
	"github.com/hhcho/sfgwas-lmm/mpc"

	"github.com/hhcho/sfgwas-lmm/crypto"
	"github.com/hhcho/sfgwas-lmm/gwas"
	"github.com/ldsec/lattigo/v2/ckks"
	"gonum.org/v1/gonum/mat"
)

type HorCipherVector []crypto.CipherVector
type HorCipherMatrix []crypto.CipherMatrix

type REGENIE struct {
	general            *gwas.ProtocolInfoLMM
	blockToChr         []int
	useADMM            bool
	B                  int
	R                  int
	Q                  int
	K                  int
	C                  int
	N                  int
	M                  int
	AllN               int
	num_iteration_lvl1 []int
	num_iteration_lvl0 []int
	lvl1_refresh_rate  int
	NProc              int

	invStd          []float64
	invStdEncrBlock []crypto.CipherVector

	Z_folds       []*mat.Dense
	Z_kminus      []*mat.Dense
	ZCache        []crypto.PlainMatrix
	ZkminusCache  []crypto.PlainMatrix
	tilde_pheno   crypto.CipherVector
	tilde_pheno_k []crypto.CipherVector //for the different folds

	ScaledZTZ          crypto.CipherMatrix
	ScaledZTZinvTwice  crypto.CipherMatrix
	ScaledZTZinvSS     mpc_core.RMat
	ScaledZTZinvSqrtSS mpc_core.RMat

	ScaledXTZ          []crypto.CipherMatrix
	ScaledXTZSS        []mpc_core.RMat
	ScaledX_iTZ_i      []*mat.Dense
	ScaledX_iTZ_iCache []crypto.PlainMatrix

	Wkarr    []crypto.CipherMatrix
	WkMean   crypto.CipherVector
	WkInvSig crypto.CipherVector

	yLOCO []HorCipherVector
}

func (reg *REGENIE) Level0() {
	m := make(map[int][]int)
	for b := 0; b < reg.B; b++ {
		folds := make([]int, reg.K)
		for k := 0; k < reg.K; k++ {
			folds[k] = k
		}
		m[b] = folds
	}
	reg.RunLevel0OnSpecificBlocks(m)
	reg.LoadAndProcessWkarr()
}
func (reg *REGENIE) GetThreadParallelMPC(offset, perBlock int) *mpc.ParallelMPC {
	mpcObjs := reg.general.GetParallelMPC()
	threadMpcObjs := make(mpc.ParallelMPC, perBlock)
	log.LLvl1("GetThreadParallelMPC", offset, perBlock, len(threadMpcObjs), len(mpcObjs))
	for p := 0; p < perBlock; p++ {
		threadMpcObjs[p] = mpcObjs[offset*perBlock+p]
	}
	return &threadMpcObjs
}
func (reg *REGENIE) LoadAndProcessWkarr() {
	cps := reg.general.GetCPS()
	mpcObj := reg.general.GetMPC()
	pid := mpcObj.GetPid()

	WkarrPar := make([][]crypto.CipherMatrix, reg.B)
	for b := 0; b < reg.B; b++ {
		WkarrPar[b] = make([]crypto.CipherMatrix, reg.K)
		for k := 0; k < reg.K; k++ {
			WkarrPar[b][k] = make(crypto.CipherMatrix, reg.R)
		}
	}

	nproc := reg.general.GetMpcMainThreads()

	jobChannels := make([]chan int, nproc)
	for i := range jobChannels {
		jobChannels[i] = make(chan int, 32)
	}

	// Dispatcher
	go func() {
		for b := 0; b < reg.B; b++ {
			jobChannels[b%nproc] <- b
		}
		for _, c := range jobChannels {
			close(c)
		}
	}()

	// Load back from disk
	var workerGroup sync.WaitGroup
	for thread := 0; thread < nproc; thread++ {
		workerGroup.Add(1)
		go func(thread int, cps *crypto.CryptoParams, mpcObj *mpc.MPC) {
			defer workerGroup.Done()
			for b := range jobChannels[thread] {
				log.LLvl1(time.Now().Format(time.StampMilli), "thread", thread, "processing block", b)
				for k := 0; k < reg.K; k++ {
					// Check if output already exists
					fname := "cipher_W_b_r" + strconv.Itoa(pid) + "_" + strconv.Itoa(k) + "_" + strconv.Itoa(b) + ".txt"
					var W_b_r crypto.CipherMatrix
					if reg.general.FileExistsForAll(mpcObj, reg.general.OutFile(fname)) {
						log.LLvl1(time.Now().Format(time.StampMilli), "Output file found, skipping computation:", fname)
						if pid > 0 {
							W_b_r = crypto.LoadCipherMatrixFromFile(cps, reg.general.OutFile(fname))
						}
					} else {
						panic("NOT ENOUGH BLOCKS DONE?")
					}
					if pid > 0 {
						for i := 0; i < reg.R; i++ {
							WkarrPar[b][k][i] = W_b_r[i]
						}
					}
				}
			}
		}(thread, cps.GetThread(thread), reg.general.GetParallelMPC()[thread])
	}
	workerGroup.Wait()

	Wkarr := make([]crypto.CipherMatrix, reg.K)
	for i := 0; i < reg.K; i++ {
		Wkarr[i] = make(crypto.CipherMatrix, reg.B*reg.R)
		ind := 0
		for b := 0; b < reg.B; b++ {
			for r := 0; r < reg.R; r++ {
				Wkarr[i][ind] = WkarrPar[b][i][r]
				ind += 1
			}
		}
	}

	if pid > 0 {
		for k := 0; k < reg.K; k++ {
			log.LLvl1(time.Now().Format(time.StampMilli), "Wk ", k)
			for p := 1; p < mpcObj.GetNParty(); p++ {
				if pid == p {
					crypto.SaveCipherMatrixToFile(cps, Wkarr[k], reg.general.OutFile("cipher_Wkarr"+strconv.Itoa(p)+"_"+strconv.Itoa(k)+".txt"))
				}
			}
		}
	}

	reg.Wkarr = Wkarr
	var scaleVariance, variance crypto.CipherVector
	log.LLvl1(time.Now().Format(time.StampMilli), "calculating variance and mean of Wkarr")
	if reg.general.FileExistsForAll(mpcObj, reg.general.OutFile("cipher_WkarrInvSig.txt")) && reg.general.FileExistsForAll(mpcObj, reg.general.OutFile("cipher_WkarrMean.txt")) {
		log.LLvl1(time.Now().Format(time.StampMilli), "Output file found, skipping computation:", "WkarrInvSig.txt and WkarrMean.txt")
		if pid > 0 {
			reg.WkInvSig = crypto.LoadCipherMatrixFromFile(cps, reg.general.OutFile("cipher_WkarrInvSig.txt"))[0]
			reg.WkMean = crypto.LoadCipherMatrixFromFile(cps, reg.general.OutFile("cipher_WkarrMean.txt"))[0]
		}
	} else {
		if pid > 0 {
			allPartiesSum, allPartiesSqScaledSum := reg.kFoldSumStats(Wkarr, reg.B*reg.R, reg.N, 1/math.Sqrt(float64(reg.AllN)-1))
			mean := crypto.CMultConst(cps, allPartiesSum, 1.0/math.Sqrt(float64(reg.AllN)), false)
			mean = crypto.CMultConst(cps, mean, 1.0/math.Sqrt(float64(reg.AllN)), false)
			mean = mpcObj.Network.CollectiveBootstrapVec(cps, mean, -1)
			reg.WkMean = mean
			crypto.SaveCipherMatrixToFile(cps, crypto.CipherMatrix{mean}, reg.general.OutFile("cipher_WkarrMean.txt"))
			variance = CalculateVariance(cps, mpcObj, allPartiesSqScaledSum, mean, reg.AllN, 1)
			scaleVariance = crypto.CMultConst(cps, variance, 1.0/float64(mpcObj.GetNParty()-1), false)
			scaleVariance = mpcObj.Network.AggregateCVec(cps, scaleVariance)

		}
		log.LLvl1(time.Now().Format(time.StampMilli), "calculating variance and mean of Wkarr done")

		invsig := CipherVectorSqrtInverse(cps, mpcObj, scaleVariance, reg.B*reg.R)
		log.LLvl1(time.Now().Format(time.StampMilli), "invsig of Wkarr done")

		if pid > 0 {
			crypto.SaveCipherMatrixToFile(cps, crypto.CipherMatrix{invsig}, reg.general.OutFile("cipher_WkarrInvSig.txt"))
		}
		reg.WkInvSig = invsig
	}

	log.LLvl1(time.Now().Format(time.StampMilli), "Level 0 done")
}
func (reg *REGENIE) MSE(preds [][]crypto.CipherVector, tilde_pheno []crypto.CipherVector) []*ckks.Ciphertext {
	cps := reg.general.GetCPS()
	network := reg.general.GetNetwork()
	MSE := make([]*ckks.Ciphertext, reg.Q)
	for q := 0; q < reg.Q; q++ {
		sumRes := CZeros(cps, cps.GetSlots())
		for k := 0; k < reg.K; k++ {
			// log.LLvl1(time.Now().Format(time.StampMilli), "lengths", len(sumRes), len(preds[k][q]), len(tilde_pheno[k]))
			res := crypto.CSub(cps, preds[k][q], tilde_pheno[k])
			resErr := CDot(cps, res, res)
			sumRes = crypto.CAdd(cps, sumRes, crypto.CipherVector{resErr})
		}
		agg := network.AggregateCVec(cps, sumRes)
		agg = network.BootstrapVecAll(cps, agg)

		agg = crypto.CMultConst(cps, agg, 1/math.Sqrt(float64(reg.AllN)), false)
		agg = crypto.CMultConst(cps, agg, 1/math.Sqrt(float64(reg.AllN)), false)

		MSE[q] = agg[0]
	}
	return MSE
}
func (reg *REGENIE) Level1() []HorCipherVector {
	cps := reg.general.GetCPS()
	mpcObj := reg.general.GetMPC()
	pid := mpcObj.GetPid()
	foldSizes := reg.general.GetGenoFoldSizes()
	h2s := compute_h2s(reg.Q)
	//artificially lower from 0.99 to 0.9
	h2s[reg.Q-1] = 0.9
	lambdas := compute_lambdas(h2s, float64(reg.B*reg.R))
	log.LLvl1(time.Now().Format(time.StampMilli), "lambdas: ", lambdas)

	allExist := true
	for k := 0; k < reg.K; k++ {
		fname := "cipher_Wkarr" + strconv.Itoa(pid) + "_" + strconv.Itoa(k) + ".txt"
		if reg.general.FileExistsForAll(mpcObj, reg.general.OutFile(fname)) {
			log.LLvl1(time.Now().Format(time.StampMilli), "Output file found :", fname)
		} else {
			allExist = false
		}
	}
	if pid > 0 && reg.Wkarr == nil && allExist {
		reg.Wkarr = make([]crypto.CipherMatrix, reg.K)
		for k := 0; k < reg.K; k++ {
			fname := "cipher_Wkarr" + strconv.Itoa(pid) + "_" + strconv.Itoa(k) + ".txt"
			log.LLvl1(time.Now().Format(time.StampMilli), "Output file found, skipping Wkarr:", fname)
			reg.Wkarr[k] = crypto.LoadCipherMatrixFromFile(cps, reg.general.OutFile(fname))
		}
	}

	fname := "cipher_WkarrInvSig.txt"
	if reg.general.FileExistsForAll(mpcObj, reg.general.OutFile(fname)) {
		log.LLvl1(time.Now().Format(time.StampMilli), "Output file found, skipping WkarrInvSig:", fname)
		if pid > 0 && reg.WkInvSig == nil {
			reg.WkInvSig = crypto.LoadCipherMatrixFromFile(cps, reg.general.OutFile(fname))[0]
			reg.WkInvSig = maskFirstPlaces(cps, reg.WkInvSig, reg.B*reg.R)
			reg.WkInvSig = mpcObj.Network.BootstrapVecAll(cps, reg.WkInvSig)
		}
	}

	fname = "cipher_WkarrMean.txt"
	if reg.general.FileExistsForAll(mpcObj, reg.general.OutFile(fname)) {
		log.LLvl1(time.Now().Format(time.StampMilli), "Output file found, skipping WkMean:", fname)
		if pid > 0 && reg.WkMean == nil {
			reg.WkMean = crypto.LoadCipherMatrixFromFile(cps, reg.general.OutFile(fname))[0]
			reg.WkMean = maskFirstPlaces(cps, reg.WkMean, reg.B*reg.R)
			reg.WkMean = mpcObj.Network.BootstrapVecAll(cps, reg.WkMean)
		}
	}
	log.LLvl1(time.Now().Format(time.StampMilli), "Wkarr loaded")

	etas := make([][]crypto.CipherVector, reg.K)
	preds := make([][]crypto.CipherVector, reg.K)
	for k := 0; k < reg.K; k++ {
		etas[k] = make([]crypto.CipherVector, reg.Q)
		preds[k] = make([]crypto.CipherVector, reg.Q)
	}
	mpcObjs := reg.GetThreadParallelMPC(0, reg.general.GetMpcPerBlock())

	XTy_sum := CZeros(cps, reg.B*reg.R)
	W_kTy_arr := make([]crypto.CipherVector, reg.K)
	if pid > 0 {
		for k := 0; k < reg.K; k++ {
			log.LLvl1(time.Now().Format(time.StampMilli), "Level1 MultLazyCMatTParallel fold", k)
			W_kTy := reg.MultLazyCMatTParallel(cps, mpcObjs, reg.Wkarr[k], crypto.CipherMatrix{reg.tilde_pheno_k[k]}, reg.WkMean, reg.WkInvSig, foldSizes[k])[0]

			// reg.Wkarr times tilde y each vector will be of size N_i
			XTy_sum = crypto.CAdd(cps, XTy_sum, W_kTy)
			W_kTy_arr[k] = W_kTy
		}
	}
	log.LLvl1(time.Now().Format(time.StampMilli), "W_kTy_arr done")
	for k := 0; k < reg.K; k++ {
		var WTy crypto.CipherVector

		if pid > 0 {
			WTy = crypto.CSub(cps, XTy_sum, W_kTy_arr[k])
		}
		log.LLvl1(time.Now().Format(time.StampMilli), "Start Ridge Regression 1")
		// Check if output already exists
		fname := "etas" + strconv.Itoa(pid) + "_" + strconv.Itoa(k) + ".txt"
		if reg.general.FileExistsForAll(mpcObj, reg.general.OutFile(fname)) {

			log.LLvl1(time.Now().Format(time.StampMilli), "Output file found, skipping computation:", fname)

			if pid > 0 {
				etas[k] = crypto.LoadCipherMatrixFromFile(cps, reg.general.OutFile("cipher_"+fname))
			}

		} else {

			log.LLvl1(time.Now().Format(time.StampMilli), "Output file not found, performing computation for:", fname)
			etas[k] = reg.RidgeRegressionLevel1(cps, mpcObjs, reg.Wkarr, reg.WkMean, reg.WkInvSig, WTy, k, lambdas)
			etas[k] = mpcObjs.GetNetworks().BootstrapMatAll(cps, etas[k])

			if pid > 0 {
				for p := 1; p < mpcObj.GetNParty(); p++ {
					if p == pid {
						crypto.SaveCipherMatrixToFile(cps, etas[k], reg.general.OutFile("cipher_etas"+strconv.Itoa(p)+"_"+strconv.Itoa(k)+".txt"))
					}
				}
			}
		}

		fname = "predslvl1" + strconv.Itoa(pid) + "_" + strconv.Itoa(k) + ".txt"
		if reg.general.FileExistsForAll(mpcObj, reg.general.OutFile(fname)) {
			log.LLvl1(time.Now().Format(time.StampMilli), "Output file found, skipping computation:", fname)
			if pid > 0 {
				preds[k] = crypto.LoadCipherMatrixFromFile(cps, reg.general.OutFile("cipher_"+fname))
			}
		} else {
			if pid > 0 {
				log.LLvl1(time.Now().Format(time.StampMilli), "Output file not found, performing computation for:", fname)
				pred := reg.MultLazyCMatParallel(cps, mpcObjs, reg.Wkarr[k], etas[k], reg.WkMean, reg.WkInvSig, foldSizes[k]) // fold size
				//need to normalize with scale again.
				for q := 0; q < reg.Q; q++ {
					pred[q] = crypto.CMultConst(cps, pred[q], 1/math.Sqrt(float64(reg.AllN)), false)
					pred[q] = crypto.CMultConst(cps, pred[q], 1/math.Sqrt(lambdas[q]), false)
				}
				pred = mpcObjs.GetNetworks().BootstrapMatAll(cps, pred)
				for p := 1; p < mpcObj.GetNParty(); p++ {
					if p == pid {
						crypto.SaveCipherMatrixToFile(cps, pred, reg.general.OutFile("cipher_predslvl1"+strconv.Itoa(p)+"_"+strconv.Itoa(k)+".txt"))
					}
				}
				preds[k] = pred
			}
		}
	}

	yLOCO := make([]HorCipherVector, reg.general.GetNumChrs())
	for i := 0; i < reg.general.GetNumChrs(); i++ {
		yLOCO[i] = make(HorCipherVector, reg.K)
	}
	minPos := 0
	if pid > 0 {
		MSE := reg.MSE(preds, reg.tilde_pheno_k)
		log.LLvl1(time.Now().Format(time.StampMilli), "MSE")
		minVal := cipherToNetworkMat(cps, mpcObj, crypto.CipherMatrix{crypto.CipherVector{MSE[0]}}, 1).At(0, 0)

		for q := 1; q < reg.Q; q++ {
			mseDecr := cipherToNetworkMat(cps, mpcObj, crypto.CipherMatrix{crypto.CipherVector{MSE[q]}}, 1).At(0, 0)
			if mseDecr < minVal {
				minVal = mseDecr
				minPos = q
			}
		}

		log.LLvl1(time.Now().Format(time.StampMilli), "finish MSE comparison minPos", minPos)
		finalPred := make(HorCipherVector, reg.K)
		for k := 0; k < reg.K; k++ {
			finalPred[k] = preds[k][minPos]
			for p := 1; p < mpcObj.GetNParty(); p++ {
				fname := reg.general.OutFile("cipher_final_pred" + strconv.Itoa(p) + "_" + strconv.Itoa(k) + ".txt")
				if pid == p {
					crypto.SaveCipherMatrixToFile(cps, crypto.CipherMatrix{finalPred[k]}, fname)
				}
			}
		}
	}

	for chr := 0; chr < reg.general.GetNumChrs(); chr++ {
		parent := reg.general.OutFile("yloco")
		err := os.MkdirAll(parent, 0700)
		if err != nil && !os.IsExist(err) {
			log.Fatal(err)
		}
		log.LLvl1("ComputeEtaMask chrom", chr)
		var mask crypto.CipherVector
		isZero := false
		if pid > 0 {
			mask, isZero = ComputeEtaMask(cps, reg.blockToChr, chr, reg.B, reg.R)
		}
		for k := 0; k < reg.K; k++ {
			fname := "yloco/cipher_yloco" + strconv.Itoa(pid) + "_" + strconv.Itoa(k) + "_" + strconv.Itoa(chr) + ".txt"
			if reg.general.FileExistsForAll(mpcObj, reg.general.OutFile(fname)) {
				log.LLvl1(time.Now().Format(time.StampMilli), "Output file found, skipping computation:", fname)
				if pid > 0 {
					reg.yLOCO[chr][k] = crypto.LoadCipherMatrixFromFile(cps, reg.general.OutFile(fname))[0]
				}

			} else if pid > 0 {
				if !isZero {
					etaLOCO := crypto.CMult(cps, mask, etas[k][minPos])

					res := reg.MultLazyCMatParallel(cps, mpcObjs, reg.Wkarr[k], crypto.CipherMatrix{etaLOCO}, reg.WkMean, reg.WkInvSig, foldSizes[k])[0] // fold size
					res = mpcObjs.GetNetworks().BootstrapVecAll(cps, res)
					res = crypto.CMultConst(cps, res, 1/math.Sqrt(float64(reg.AllN)), false)
					res = crypto.CMultConst(cps, res, 1/math.Sqrt(float64(lambdas[minPos])), false)
					yLOCO[chr][k] = mpcObjs.GetNetworks().BootstrapVecAll(cps, res)
				} else {
					yLOCO[chr][k] = CZeros(cps, foldSizes[k])
				}
				for p := 1; p < mpcObj.GetNParty(); p++ {
					fname := reg.general.OutFile("yloco/cipher_yloco" + strconv.Itoa(p) + "_" + strconv.Itoa(k) + "_" + strconv.Itoa(chr) + ".txt")
					if pid == p {
						crypto.SaveCipherMatrixToFile(cps, crypto.CipherMatrix{yLOCO[chr][k]}, fname)
						log.Info(">>>>> Saved yloco: ", fname)
					}
				}
			}
		}

	}

	reg.yLOCO = yLOCO
	return reg.yLOCO
}
func (reg *REGENIE) Load_Completed_YLOCO() []HorCipherVector {
	cps := reg.general.GetCPS()
	mpcObj := reg.general.GetMPC()
	pid := mpcObj.GetPid()
	yLOCO := make([]HorCipherVector, reg.general.GetNumChrs())
	for i := 0; i < reg.general.GetNumChrs(); i++ {
		yLOCO[i] = make(HorCipherVector, reg.K)
	}
	for chr := 0; chr < reg.general.GetNumChrs(); chr++ {
		for k := 0; k < reg.K; k++ {
			fname := "yloco/cipher_yloco" + strconv.Itoa(pid) + "_" + strconv.Itoa(k) + "_" + strconv.Itoa(chr) + ".txt" // fname := "yloco/yloco" + strconv.Itoa(pid) + "_" + strconv.Itoa(k) + "_" + strconv.Itoa(chr) + ".txt"
			if reg.general.FileExistsForAll(mpcObj, reg.general.OutFile(fname)) {
				log.LLvl1(time.Now().Format(time.StampMilli), "Output file found, skipping computation:", fname)
				if pid > 0 {
					yLOCO[chr][k] = crypto.LoadCipherMatrixFromFile(cps, reg.general.OutFile(fname))[0]
				}

			} else {
				panic(fmt.Sprintf("Didn't have ylocos: %s", reg.general.OutFile(fname)))
			}
		}

	}
	reg.yLOCO = yLOCO
	return reg.yLOCO
}

// returns one CipherVector per block
func (reg *REGENIE) AssociationStats(y_hat_LOCO []HorCipherVector) []crypto.CipherVector {
	debug := true
	parent := reg.general.OutFile("assoc")
	err := os.MkdirAll(parent, 0700)
	if err != nil && !os.IsExist(err) {
		log.Fatal(err)
	}
	numChr := len(y_hat_LOCO)

	cps := reg.general.GetCPS()
	mpcObj := reg.general.GetMPC()
	pid := mpcObj.GetPid()

	yres_LOCO := make([][]crypto.CipherVector, numChr)
	for chr := range yres_LOCO {
		yres_LOCO[chr] = make([]crypto.CipherVector, reg.K)
	}

	nproc := reg.general.GetMpcMainThreads()
	nprocPerBlock := reg.general.GetLocalThreads() / nproc
	reg.NProc = nprocPerBlock
	log.LLvl1(time.Now().Format(time.StampMilli), "Using", nproc, "threads in parallel with", nprocPerBlock, "inner threads each")

	sigmaESS := mpc_core.InitRVec(mpcObj.GetRType().Zero(), numChr)
	if pid > 0 {
		// Dispatcher
		jobChannels := make([]chan int, nproc)
		for i := range jobChannels {
			jobChannels[i] = make(chan int, 32)
		}
		go func() {
			for chr := 0; chr < numChr; chr++ {
				jobChannels[chr%nproc] <- chr
			}
			for _, c := range jobChannels {
				close(c)
			}
		}()

		// Workers
		var workerGroup sync.WaitGroup
		for thread := 0; thread < nproc; thread++ {
			workerGroup.Add(1)
			go func(thread int, cps *crypto.CryptoParams, mpcObj *mpc.MPC) {
				defer workerGroup.Done()

				for chr := range jobChannels[thread] {
					log.LLvl1(time.Now().Format(time.StampMilli), "Computing yres_LOCO thread", thread, "chr", chr)

					var sum *ckks.Ciphertext
					for fold := range y_hat_LOCO[chr] {
						yres_LOCO[chr][fold] = crypto.CSub(cps, reg.tilde_pheno_k[fold], y_hat_LOCO[chr][fold])
						yres_LOCO_sq := crypto.CMult(cps, yres_LOCO[chr][fold], yres_LOCO[chr][fold])
						res := crypto.InnerSumAll(cps, yres_LOCO_sq)
						if fold == 0 {
							sum = res
						} else {
							sum = crypto.Add(cps, sum, res)
						}
					}

					sum = mpcObj.Network.AggregateCText(cps, sum)

					inter := crypto.CMultConst(cps, crypto.CipherVector{sum}, 1/math.Sqrt(float64(reg.AllN-reg.C+1)), false)
					sigmaE := crypto.CMultConst(cps, inter, 1/math.Sqrt(float64(reg.AllN-reg.C+1)), false)[0]

					sigmaESS[chr] = mpcObj.CiphertextToSS(cps, mpcObj.GetRType(), sigmaE, -1, 1)[0]
					yres_LOCO[chr] = mpcObj.Network.BootstrapMatAll(cps, yres_LOCO[chr])
				}

			}(thread, cps.GetThread(thread), reg.general.GetParallelMPC()[thread])
		}
		workerGroup.Wait()
	}

	_, sigmaEInvSS := mpcObj.SqrtAndSqrtInverse(sigmaESS)
	sigmaEInv := make(crypto.CipherVector, numChr)
	if pid > 0 {
		for chr := 0; chr < numChr; chr++ {
			sigmaEInv[chr] = mpcObj.SStoCiphertext(cps, mpc_core.RVec{sigmaEInvSS[chr]})
			sigmaEInv[chr] = crypto.Rebalance(cps, sigmaEInv[chr])
		}
	}

	var RT mpc_core.RMat
	if pid > 0 {
		RT = reg.ScaledZTZinvSqrtSS
	} else {
		RT = mpc_core.InitRMat(mpcObj.GetRType().Zero(), reg.C, reg.C)
	}

	// Dispatcher
	jobChannels := make([]chan int, nproc)
	for i := range jobChannels {
		jobChannels[i] = make(chan int, 32)
	}
	go func() {
		for chr := 0; chr < numChr; chr++ {
			jobChannels[chr%nproc] <- chr
		}
		for _, c := range jobChannels {
			close(c)
		}
	}()

	w := make([][]crypto.CipherVector, numChr)

	// Workers
	var workerGroup sync.WaitGroup
	for thread := 0; thread < nproc; thread++ {
		workerGroup.Add(1)
		go func(thread int, cps *crypto.CryptoParams, mpcObj *mpc.MPC) {
			defer workerGroup.Done()

			for chr := range jobChannels[thread] {
				log.LLvl1(time.Now().Format(time.StampMilli), "Computing Py_LOCO thread", thread, "chr", chr)
				w[chr] = make([]crypto.CipherVector, reg.K)

				var ZTMul crypto.CipherVector
				if pid > 0 {
					var Z_iTMul crypto.CipherVector
					for fold := 0; fold < reg.K; fold++ {
						res := reg.MulZTCache(cps, crypto.CipherMatrix{yres_LOCO[chr][fold]}, fold)[0]
						if fold == 0 {
							Z_iTMul = res
						} else {
							Z_iTMul = crypto.CAdd(cps, Z_iTMul, res)
						}
					}
					ZTMul = mpcObj.Network.AggregateCVec(cps, Z_iTMul)
				}

				ZTMulSS := mpcObj.CVecToSS(cps, mpcObj.GetRType(), ZTMul, -1, 1, reg.C)

				var ScaledZTZinvSS mpc_core.RMat
				if pid > 0 {
					ScaledZTZinvSS = reg.ScaledZTZinvSS
				} else {
					ScaledZTZinvSS = mpc_core.InitRMat(mpcObj.GetRType().Zero(), reg.C, reg.C)
				}

				ZTZinv_ZTMulSS := mpcObj.SSMultMat(mpc_core.RMat{ZTMulSS}, ScaledZTZinvSS)
				ZTZinv_ZTMulSS = mpcObj.TruncMat(ZTZinv_ZTMulSS, mpcObj.GetDataBits(), mpcObj.GetFracBits())

				ZTZinv_ZTMul := crypto.CipherMatrix{mpcObj.SSToCVec(cps, ZTZinv_ZTMulSS[0])}

				if pid > 0 {
					ZTZinv_ZTMul[0][0] = crypto.MaskTrunc(cps, ZTZinv_ZTMul[0][0], reg.C)
				}
				if pid > 0 {
					Z_ZTZinv_ZTMul := make([]crypto.CipherVector, reg.K)
					for fold := 0; fold < reg.K; fold++ {
						Z_ZTZinv_ZTMul[fold] = reg.MulZCache(cps, ZTZinv_ZTMul, fold, false)[0]
						tmp := crypto.CMultConst(cps, Z_ZTZinv_ZTMul[fold], 1/math.Sqrt(float64(reg.AllN)), false)
						w[chr][fold] = crypto.CSub(cps, yres_LOCO[chr][fold], tmp)
					}
				}
			}
		}(thread, cps.GetThread(thread), reg.general.GetParallelMPC()[thread])
	}
	workerGroup.Wait()

	/* Compute association stats */

	// Dispatcher
	for i := range jobChannels {
		jobChannels[i] = make(chan int, 32)
	}
	go func() {
		for block := 0; block < reg.B; block++ {
			jobChannels[block%nproc] <- block
		}
		for _, c := range jobChannels {
			close(c)
		}
	}()

	stats := make([]crypto.CipherVector, reg.B)

	// Workers
	for thread := 0; thread < nproc; thread++ {
		workerGroup.Add(1)
		go func(thread int, cps *crypto.CryptoParams, mpcObj *mpc.MPC) {
			defer workerGroup.Done()

			for block := range jobChannels[thread] {
				log.LLvl1(time.Now().Format(time.StampMilli), "Computing final statistics, thread", thread, "block", block)

				// Compute g'g (sqSum), Z'g, and w'g
				var wTgLocal crypto.CipherMatrix
				var gTgPlain []float64

				if pid > 0 {
					gTgPlain = make([]float64, reg.general.GetGenoBlockSizes()[block])

					for fold := 0; fold < reg.K; fold++ {
						log.LLvl1(time.Now().Format(time.StampMilli), "w", fold, block, w[reg.blockToChr[block]][fold])

						res, sqSum := reg.MulXTStreamFold(cps, crypto.CipherMatrix{w[reg.blockToChr[block]][fold]}, fold, block, true)

						log.LLvl1(time.Now().Format(time.StampMilli), "MulXTStream fold", fold, "len(sqSum)", len(sqSum), "finished")

						if fold == 0 {
							wTgLocal = res
						} else {
							wTgLocal = CMatAdd(cps, wTgLocal, res)
						}

						if sqSum != nil {
							if len(sqSum) != len(gTgPlain) {
								panic("len(sqSum) != number of variants in block")
							}

							for j := range gTgPlain {
								gTgPlain[j] += sqSum[j]
							}
						}
					}
				}

				// Aggregate across parties
				var gPy crypto.CipherVector
				if pid > 0 {
					gPy = mpcObj.Network.AggregateCVec(cps, wTgLocal[0])
				}
				gPy = mpcObj.Network.CollectiveBootstrapVec(cps, gPy, -1)

				var zTg mpc_core.RMat
				if pid > 0 {
					zTg = reg.ScaledXTZSS[block] // precomputed
				} else {
					zTg = mpc_core.InitRMat(mpcObj.GetRType().Zero(), reg.C, reg.general.GetGenoBlockSizes()[block])
				}

				// Compute u = R'Z'g
				log.LLvl1(time.Now().Format(time.StampMilli), "RT dims", len(RT), len(RT[0]))
				log.LLvl1(time.Now().Format(time.StampMilli), "zTg dims", len(zTg), len(zTg[0]))

				u := mpcObj.SSMultMat(RT, zTg)
				u = mpcObj.TruncMat(u, mpcObj.GetDataBits(), mpcObj.GetFracBits())

				// Compute u'u
				uSq := mpcObj.SSMultElemMat(u, u)
				uSqSum := uSq.Sum(1)
				uSqSum = mpcObj.TruncVec(uSqSum, mpcObj.GetDataBits(), mpcObj.GetFracBits())
				if debug && block == 0 {
					uSqSumPlain := mpcObj.RevealSymVec(uSqSum).ToFloat(mpcObj.GetFracBits())
					gwas.SaveFloatVectorToFile(reg.general.OutFile("assoc/uSqSumPlain.txt"), uSqSumPlain)
					gwas.SaveFloatVectorToFile(reg.general.OutFile("assoc/gTgPlain.txt"), gTgPlain)
				}
				// Compute g'Pg = g'g - u'u
				var gPg mpc_core.RVec
				if pid > 0 {
					gPg = mpc_core.FloatToRVec(mpcObj.GetRType(), gTgPlain, mpcObj.GetFracBits())

					log.LLvl1(time.Now().Format(time.StampMilli), "gPg dims", len(gPg))
					log.LLvl1(time.Now().Format(time.StampMilli), "uSqSum dims", len(uSqSum))

					gPg.Sub(uSqSum)
				} else {
					gPg = mpc_core.InitRVec(mpcObj.GetRType().Zero(), reg.general.GetGenoBlockSizes()[block])
				}

				_, gPgInvSqrtSS := mpcObj.SqrtAndSqrtInverse(gPg)

				gPgInvSqrt := mpcObj.SSToCVec(cps, gPgInvSqrtSS)

				if debug && block == 0 {
					gPgPlain := mpcObj.RevealSymVec(gPg).ToFloat(mpcObj.GetFracBits())
					gPgInvSqrtSSPlain := mpcObj.RevealSymVec(gPgInvSqrtSS).ToFloat(mpcObj.GetFracBits())

					gwas.SaveFloatVectorToFile(reg.general.OutFile("assoc/gPg.txt"), gPgPlain)
					gwas.SaveFloatVectorToFile(reg.general.OutFile("assoc/gPgInvSqrtSS.txt"), gPgInvSqrtSSPlain)
				}

				if pid > 0 {
					stats[block] = crypto.CMult(cps, gPy, gPgInvSqrt)
					stats[block] = crypto.CMultScalar(cps, stats[block], sigmaEInv[reg.blockToChr[block]])
				}

			}
		}(thread, cps.GetThread(thread), reg.general.GetParallelMPC()[thread])
	}
	workerGroup.Wait()

	return stats
}
func (reg *REGENIE) generateBlockCache(fname string, p, k, b int) string {
	return reg.general.CacheFile(gwas.GenerateFName(fname, p, k, b))
}
func (reg *REGENIE) SaveMatDenseToFile(x *mat.Dense, fname string, k, b int) {
	mpcObj := reg.general.GetMPC()
	parent := reg.general.CacheFile("block_" + strconv.Itoa(b))
	err := os.MkdirAll(parent, 0700)
	if err != nil && !os.IsExist(err) {
		log.Fatal(err)
	}
	for p := 1; p < mpcObj.GetNParty(); p++ {
		gwas.SaveMatDenseToFile(mpcObj, p, x, reg.generateBlockCache(fname, p, k, b))
	}
}
func (reg *REGENIE) GetBlockSizeInfo(b int) (int, int) {
	mpcObj := reg.general.GetMPC()
	pid := mpcObj.GetPid()
	blockSize := 0
	if pid == 0 {
		return 0, 0
	}
	blockSizes := reg.general.GetGenoBlockSizes()
	blockSize = blockSizes[b]
	runningSum := 0
	for i := 0; i < b; i++ {
		runningSum += blockSizes[i]
	}
	return blockSize, runningSum
}
func (reg *REGENIE) RunLevel0OnSpecificBlocks(spec map[int][]int) {
	reg.NProc = (reg.general.GetLocalThreads() / reg.general.GetMpcMainThreads())
	spec_blocks := make([]int, 0)
	for k := range spec {
		spec_blocks = append(spec_blocks, k)
	}
	sort.Ints(spec_blocks)
	parNet := reg.general.GetParallelNetwork()
	nproc := reg.general.GetMpcMainThreads()
	firstStart := time.Now()
	cps := reg.general.GetCPS()
	mpcObj := reg.general.GetMPC()
	pid := mpcObj.GetPid()
	blockSizes := reg.general.GetGenoBlockSizes()
	h2s := compute_h2s(reg.R)
	lambdas := compute_lambdas(h2s, float64(reg.M))
	log.LLvl1(time.Now().Format(time.StampMilli), "lambdas: ", lambdas)

	jobChannels := make([]chan int, nproc)
	for i := range jobChannels {
		jobChannels[i] = make(chan int, 32)
	}
	// Dispatcher
	go func() {
		for spec_idx := 0; spec_idx < len(spec_blocks); spec_idx++ { // main difference
			b := spec_blocks[spec_idx]
			jobChannels[b%nproc] <- b
		}
		for _, c := range jobChannels {
			close(c)
		}
	}()
	calcGTicketsChannel := make(chan int, reg.general.GetCalcGCapacity())
	for i := 0; i < reg.general.GetCalcGCapacity(); i++ { // Make space for concurrent runs
		calcGTicketsChannel <- 1
	}

	// Workers
	var workerGroup sync.WaitGroup
	for thread := 0; thread < nproc; thread++ {
		workerGroup.Add(1)
		go func(thread int, cps *crypto.CryptoParams, mpcObjs *mpc.ParallelMPC) {
			mpcObj := (*mpcObjs)[0]
			defer workerGroup.Done()
			for b := range jobChannels[thread] {
				log.LLvl1(time.Now().Format(time.StampMilli), "thread", thread, "processing block", b)
				spec_folds := spec[b]
				for _, k := range spec_folds {
					log.LLvl1(time.Now().Format(time.StampMilli), "thread", thread, "processing fold", k)
					W_b_r_fname := reg.general.OutFile("W_b_r" + strconv.Itoa(pid) + "_" + strconv.Itoa(k) + "_" + strconv.Itoa(b) + ".txt")
					var W_b_r crypto.CipherMatrix
					if reg.general.FileExistsForAll(mpcObj, W_b_r_fname) {
						log.LLvl1(time.Now().Format(time.StampMilli), "Output file found, skipping computation:", W_b_r_fname)
					} else {
						ATAs := make([]*mat.Dense, reg.K)
						ScaledGenoTZs := make([]*mat.Dense, reg.K)
						XTys := make([]crypto.CipherVector, reg.K)
						var ATASum *mat.Dense
						var ScaledGenoTZSum *mat.Dense
						var XTySum crypto.CipherVector
						if pid > 0 && reg.useADMM {
							ATASum = mat.NewDense(blockSizes[b], blockSizes[b], nil)
							ScaledGenoTZSum = mat.NewDense(blockSizes[b], reg.C, nil)
							XTySum = CZeros(cps, blockSizes[b])
							//compute all ATA for each fold keeping block b constant
							log.LLvl1(time.Now().Format(time.StampMilli), "Block", b, "Calculating ATA overhead")
							for preK := 0; preK < reg.K; preK++ {
								ATAfname := reg.generateBlockCache("ATA", pid, preK, b)
								ScaledGenoTZfname := reg.generateBlockCache("ScaledGenoTZs", pid, preK, b)
								XTyfname := reg.generateBlockCache("XTy", pid, preK, b)
								if reg.general.Exists(ATAfname) && reg.general.Exists(ScaledGenoTZfname) && reg.general.Exists(XTyfname) {
									ATAs[preK] = gwas.LoadMatDenseCacheFromFile(ATAfname)
									r, c := ATAs[preK].Dims()
									if r != blockSizes[b] || c != blockSizes[b] {
										log.Fatal("Matrix dims not expected")
									}
									ScaledGenoTZs[preK] = gwas.LoadMatDenseCacheFromFile(ScaledGenoTZfname)
									r, c = ScaledGenoTZs[preK].Dims()
									if r != blockSizes[b] || c != reg.C {
										log.Fatal("Matrix dims not expected")
									}
									XTys[preK] = gwas.LoadCacheFromFile(cps, XTyfname)[0]

									ScaledGenoTZSum.Add(ScaledGenoTZSum, ScaledGenoTZs[preK])
									ATASum.Add(ATASum, ATAs[preK])
									XTySum = crypto.CAdd(cps, XTySum, XTys[preK])

									log.LLvl1(time.Now().Format(time.StampMilli), "Block", b, preK, "Found ATA and ScaledGenoTZSum")
									continue
								}
								blockSize, runningSum := reg.GetBlockSizeInfo(b)
								//calculate invstd for this region
								invStdSlice := make([]float64, blockSize)
								copy(invStdSlice, reg.invStd[runningSum:runningSum+blockSize])
								log.LLvl1(time.Now().Format(time.StampMilli), "Block", b, "ATA big mult")
								ATA := gwas.MatSATASClear(reg.general.GetGeno()[preK][b], invStdSlice)
								log.LLvl1(time.Now().Format(time.StampMilli), "Block", b, "ATA after big mult")
								ATASum.Add(ATASum, ATA)
								ATAs[preK] = ATA

								ScaledGenoTZ := gwas.MatTSClearSimple(reg.general.GetGeno()[preK][b], reg.Z_folds[preK], invStdSlice)
								ScaledGenoTZ.Scale(1/math.Sqrt(float64(reg.AllN)), ScaledGenoTZ)
								ScaledGenoTZSum.Add(ScaledGenoTZSum, ScaledGenoTZ)
								ScaledGenoTZs[preK] = ScaledGenoTZ

								pheno := crypto.CopyEncryptedVector(reg.tilde_pheno_k[preK])
								XTy := reg.XtildeTMulVecStreamLocal(cps, mpcObj, pheno, preK, b)
								XTySum = crypto.CAdd(cps, XTySum, XTy)
								XTys[preK] = XTy

								reg.SaveMatDenseToFile(ATAs[preK], "ATA", preK, b)
								reg.SaveMatDenseToFile(ScaledGenoTZs[preK], "ScaledGenoTZs", preK, b)
								for p := 1; p < mpcObj.GetNParty(); p++ {
									fname := reg.generateBlockCache("XTy", p, preK, b)
									gwas.SaveMatrixToFile(cps, mpcObj, crypto.CipherMatrix{XTy}, blockSize, p, fname)
								}
							}
							log.LLvl1(time.Now().Format(time.StampMilli), "Block", b, "ATA overhead done")
						}
						var Z_iTZ_iZTZinv crypto.CipherMatrix
						if pid > 0 && reg.useADMM {
							log.LLvl1(time.Now().Format(time.StampMilli), "Block", b, "Fold", k, "Calculating Z_iTZ_iZTZinv overhead")
							//precompute Z and X components
							Z_kminus := reg.Z_kminus[k]
							_, c := Z_kminus.Dims()
							ScaledZ_kminusTZ_kminus := mat.NewDense(c, c, nil)
							ScaledZ_kminusTZ_kminus.Mul(Z_kminus.T(), Z_kminus)
							ScaledZ_kminusTZ_kminusPlain := crypto.EncodeDense(cps, ScaledZ_kminusTZ_kminus)
							EvenBiggerInv := reg.ScaledZTZinvTwice
							Z_iTZ_iZTZinv = gwas.CPMatMult5(cps, ScaledZ_kminusTZ_kminusPlain, EvenBiggerInv) //should have no scale
							Z_iTZ_iZTZinv = mpcObj.Network.BootstrapMatAll(cps, Z_iTZ_iZTZinv)
							Z_iTZ_iZTZinv = CMatMultConst(cps, Z_iTZ_iZTZinv, 1/math.Sqrt(float64(reg.AllN)))
							Z_iTZ_iZTZinv = mpcObj.Network.BootstrapMatAll(cps, Z_iTZ_iZTZinv)
						}

						log.LLvl1(time.Now().Format(time.StampMilli), "Block", b, "Fold", k, "Calculating predictors")
						start := time.Now()

						var predictors crypto.CipherMatrix
						if reg.useADMM {
							var ATA_k *mat.Dense
							var GenoTZ_k *mat.Dense
							var XTy_k crypto.CipherVector
							if pid > 0 {
								ATA_k = mat.NewDense(blockSizes[b], blockSizes[b], nil)
								ATA_k.Sub(ATASum, ATAs[k])
								GenoTZ_k = mat.NewDense(blockSizes[b], reg.C, nil)
								GenoTZ_k.Sub(ScaledGenoTZSum, ScaledGenoTZs[k])
								XTy_k = crypto.CSub(cps, XTySum, XTys[k])
							}
							predictors = reg.RidgeRegressionADMMWood(cps, mpcObjs, lambdas, Z_iTZ_iZTZinv, ATA_k, GenoTZ_k, XTy_k, k, b, calcGTicketsChannel)
						} else {
							predictors = reg.RidgeRegressionCGD(cps, mpcObj, lambdas, k, b)
						}

						log.LLvl1(time.Now().Format(time.StampMilli), "Block", b, "Fold", k, "Calculating predictors done", time.Since(start))

						if pid > 0 {
							log.LLvl1(time.Now().Format(time.StampMilli), "Block", b, "Fold", k, "Saving predictors")
							for p := 1; p < mpcObj.GetNParty(); p++ {
								if p == pid {
									crypto.SaveCipherMatrixToFile(cps, predictors, reg.general.OutFile("cipher_predictors"+strconv.Itoa(p)+"_"+strconv.Itoa(k)+"_"+strconv.Itoa(b)+".txt"))
								}
							}
							//bottom step is expensive
							W_b_r = reg.XtildeMulStream(cps, mpcObj, predictors, k, b, false)
							for p := 1; p < mpcObj.GetNParty(); p++ {
								if p == pid {
									crypto.SaveCipherMatrixToFile(cps, W_b_r, reg.general.OutFile("cipher_W_b_r"+strconv.Itoa(p)+"_"+strconv.Itoa(k)+"_"+strconv.Itoa(b)+".txt"))
								}
							}
							log.LLvl1(time.Now().Format(time.StampMilli), "Block", b, "Fold", k, "Saved W_b_r")
						}
					}
				}
			}
		}(thread, cps.GetThread(thread), reg.GetThreadParallelMPC(thread, reg.general.GetMpcPerBlock()))
	}
	workerGroup.Wait()

	log.LLvl1(time.Now().Format(time.StampMilli), "Level 0 subset of blocks done")
	log.LLvl1(time.Now().Format(time.StampMilli), "Level0 time: ", time.Now().Sub(firstStart).String())
	log.LLvl1(time.Now().Format(time.StampMilli), "Level0 Network Log")
	parNet.PrintNetworkLog()
	parNet.ResetNetworkLog()
}
func (reg *REGENIE) RunREGENIE() {
	yLOCO := reg.RunStep1()
	reg.RunStep2(yLOCO)
}
func (reg *REGENIE) RunStep1() []HorCipherVector {
	parNet := reg.general.GetParallelNetwork()
	mpcObj := reg.general.GetMPC()
	start := time.Now()
	reg.PrecomputeVals(true)
	log.LLvl1(time.Now().Format(time.StampMilli), "Preprocessing time: ", time.Now().Sub(start).String())

	log.LLvl1(time.Now().Format(time.StampMilli), "Precompute Network Log")
	parNet.PrintNetworkLog()
	parNet.ResetNetworkLog()

	start = time.Now()
	reg.Level0()
	log.LLvl1(time.Now().Format(time.StampMilli), "Level0 time: ", time.Now().Sub(start).String())

	log.LLvl1(time.Now().Format(time.StampMilli), "Level0 Network Log")
	parNet.PrintNetworkLog()
	parNet.ResetNetworkLog()

	start = time.Now()
	yLOCO := reg.Level1()
	log.LLvl1(time.Now().Format(time.StampMilli), "Level1 time: ", time.Now().Sub(start).String())

	log.LLvl1(time.Now().Format(time.StampMilli), "Level1 Network Log")
	parNet.PrintNetworkLog()
	parNet.ResetNetworkLog()
	mpcObj.AssertSync()
	return yLOCO
}
func (reg *REGENIE) RunStep2(yLOCO []HorCipherVector) {
	parNet := reg.general.GetParallelNetwork()
	mpcObj := reg.general.GetMPC()
	reg.PrecomputeVals(false)

	start := time.Now()

	reg.AssociationStats(yLOCO)
	log.LLvl1(time.Now().Format(time.StampMilli), "AssociationStats time: ", time.Now().Sub(start).String())

	log.LLvl1(time.Now().Format(time.StampMilli), "AssociationStats Network Log")
	parNet.PrintNetworkLog()
	parNet.ResetNetworkLog()

	mpcObj.AssertSync()
}
func (reg *REGENIE) kFoldSumStats(X []crypto.CipherMatrix, vecLen int, n int, scale float64) (crypto.CipherVector, crypto.CipherVector) {
	network := reg.general.GetNetwork()
	cps := reg.general.GetCPS()
	allSum := CZeros(cps, vecLen)
	allSqScaledSum := CZeros(cps, vecLen)
	for k := 0; k < reg.K; k++ {
		sum, sqScaledSum := sumStats(cps, X[k], n, scale)
		allSum = crypto.CAdd(cps, allSum, sum)
		allSqScaledSum = crypto.CAdd(cps, allSqScaledSum, sqScaledSum)
	}
	allPartiesSum := network.AggregateCVec(cps, allSum)
	allPartiesSqScaledSum := network.AggregateCVec(cps, allSqScaledSum)
	return allPartiesSum, allPartiesSqScaledSum
}
func (reg *REGENIE) ComputeStdInv() *mat.Dense {
	if reg.general.GetMPC().GetPid() == 0 {
		return nil
	}
	log.LLvl1("reg.M", reg.M)
	pid := reg.general.GetMPC().GetPid()
	cps := reg.general.GetCPS()
	mpcObj := reg.general.GetMPC()

	fname := "stdInv_" + strconv.Itoa(pid) + ".txt"
	if reg.general.CacheExists(fname) {
		log.LLvl1("Stdev cache found")
	} else {
		log.LLvl1("Calculating stdev")

		// Compute local sum(x) and sum(x^2)
		sx := make([]float64, reg.M)
		sx2 := make([]float64, reg.M)
		for f := 0; f < reg.K; f++ {
			nind := reg.general.GetGenoFoldSizes()[f]
			shift := 0
			for b := 0; b < reg.B; b++ {
				nsnp := reg.general.GetGenoBlockSizes()[b]

				gfs := reg.general.GetGeno()[f][b]
				gfs.Reset()
				for r := 0; r < nind; r++ {
					row := gfs.NextRow()
					if len(row) != nsnp {
						panic(fmt.Sprint("GenoFileStream has inconsistent number of snps: fold", f, "block", b, "expected", nsnp, "got", len(row)))
					}

					for c := range row {
						sx[shift+c] += float64(row[c])
						sx2[shift+c] += float64(row[c] * row[c])
					}
				}

				shift += nsnp
			}
		}

		sxEnc, _ := crypto.EncryptFloatVector(cps, sx)
		sx2Enc, _ := crypto.EncryptFloatVector(cps, sx2)

		sxEnc = mpcObj.Network.AggregateCVec(cps, sxEnc)
		sx2Enc = mpcObj.Network.AggregateCVec(cps, sx2Enc)

		sxEnc = crypto.CMultConst(cps, sxEnc, 1/math.Sqrt(float64(reg.AllN-reg.C+1)), true) // sx / sqrt(N)
		sxEnc = crypto.CMult(cps, sxEnc, sxEnc)                                             // (sx)^2 / N

		varEnc := crypto.CSub(cps, sx2Enc, sxEnc)
		varEnc = crypto.CMultConst(cps, varEnc, 1/float64(reg.AllN-reg.C+1), true)

		varPlain := mpcObj.Network.CollectiveDecryptVec(cps, varEnc, -1)

		varFloat := crypto.DecodeFloatVector(cps, varPlain)[:reg.M]
		varInvSqrtFloat := make([]float64, len(varFloat))
		for i := 0; i < len(varFloat); i++ {
			varInvSqrtFloat[i] = 1 / math.Sqrt(varFloat[i])
		}

		log.LLvl1("First few inv stdev:", varInvSqrtFloat[:10])
		log.LLvl1("reg.M", reg.M)

		gwas.SaveFloatVectorToFile(reg.general.CacheFile(fname), varInvSqrtFloat)

		log.LLvl1("Stdev saved to cache")

		varInvSqrtFloatDense := mat.NewDense(len(varInvSqrtFloat), 1, varInvSqrtFloat)
		return varInvSqrtFloatDense
	}

	return gwas.LoadMatrixFromFile(reg.general.CacheFile(fname), '\t')
}
func (reg *REGENIE) LoadGFS(isStep1 bool) {
	config := reg.general.GetConfig()
	mpcObj := reg.general.GetMPC()
	pid := mpcObj.GetPid()

	reg.K = config.GenoNumFolds
	reg.C = config.NumCovs
	reg.N = config.NumInds[pid]
	totalInds := 0
	for i := range config.NumInds {
		totalInds += config.NumInds[i]
	}
	reg.AllN = totalInds
	genoFoldSizes := readIntSliceFromFile(config.GenoFoldSizeFile, "geno_fold_size_file", pid, config.GenoNumFolds)
	reg.general.SetGenoFoldSizes(genoFoldSizes)

	prefix := ""
	if isStep1 {
		reg.B = config.GenoNumBlocks
		reg.M = config.NumSnps
		genoBlockSizes := readIntSliceFromFile(config.GenoBlockSizeFile, "geno_block_size_file", pid, config.GenoNumBlocks)
		genoBlockToChr := readIntSliceFromFile(config.GenoBlockToChromFile, "geno_block_to_chrom_file", pid, config.GenoNumBlocks)
		numChrs := maxIntSlice(genoBlockToChr) + 1

		if sumIntSlice(genoBlockSizes) != config.NumSnps {
			log.Fatalf("Sum of block sizes does not match number of snps")
		}

		if pid > 0 {
			if sumIntSlice(genoFoldSizes) != config.NumInds[pid] {
				log.Fatalf("Sum of fold sizes does not match number of inds")
			}
		}
		reg.blockToChr = genoBlockToChr
		reg.general.SetGenoBlockToChr(genoBlockToChr)
		reg.general.SetGenoBlockSizes(genoBlockSizes)
		reg.general.SetNumChrs(numChrs)

		gwasParams := gwas.InitGWASParams(config.NumInds, config.NumSnps, config.NumCovs, config.NumPCs, config.SnpDistThres)
		reg.general.SetGWASParams(gwasParams)
		if pid > 0 {
			pos := gwas.LoadSNPPositionFile(config.SnpPosFile, '\t')
			log.LLvl1(time.Now(), "First few SNP positions:", pos[:5])
			reg.general.SetPos(pos)
		}
		prefix = reg.general.GetConfig().GenoBinFilePrefix

	} else {
		reg.B = config.Step2GenoNumBlocks
		reg.M = config.Step2NumSnps
		genoBlockSizes := readIntSliceFromFile(config.Step2GenoBlockSizeFile, "step_2_geno_block_size_file", pid, config.Step2GenoNumBlocks)
		genoBlockToChr := readIntSliceFromFile(config.Step2GenoBlockToChromFile, "step_2_geno_block_to_chrom_file", pid, config.Step2GenoNumBlocks)
		numChrs := maxIntSlice(genoBlockToChr) + 1

		if sumIntSlice(genoBlockSizes) != config.NumSnps {
			log.Fatalf("Sum of block sizes does not match number of snps")
		}

		if pid > 0 {
			if sumIntSlice(genoFoldSizes) != config.NumInds[pid] {
				log.Fatalf("Sum of fold sizes does not match number of inds")
			}
		}
		reg.blockToChr = genoBlockToChr
		reg.general.SetGenoBlockToChr(genoBlockToChr)
		reg.general.SetGenoBlockSizes(genoBlockSizes)
		reg.general.SetNumChrs(numChrs)

		gwasParams := gwas.InitGWASParams(config.NumInds, config.Step2NumSnps, config.NumCovs, config.NumPCs, config.SnpDistThres)
		reg.general.SetGWASParams(gwasParams)
		if pid > 0 {
			pos := gwas.LoadSNPPositionFile(config.Step2SnpPosFile, '\t')
			log.LLvl1(time.Now(), "First few SNP positions:", pos[:5])
			reg.general.SetPos(pos)
		}
		prefix = reg.general.GetConfig().Step2GenoBinFilePrefix
	}

	if pid > 0 {
		genofs := make([][]*gwas.GenoFileStream, reg.K)
		genoTfs := make([][]*gwas.GenoFileStream, reg.K)
		for i := 0; i < reg.K; i++ {
			genofs[i] = make([]*gwas.GenoFileStream, reg.B)
			genoTfs[i] = make([]*gwas.GenoFileStream, reg.B)
		}
		var blockSizes, foldSizes []int
		if pid > 0 {
			blockSizes = reg.general.GetGenoBlockSizes()
			foldSizes = reg.general.GetGenoFoldSizes()
		}
		for i := 0; i < reg.K; i++ {
			for j := 0; j < reg.B; j++ {
				rows := foldSizes[i]
				foldSizes[i] = rows
				filename := fmt.Sprintf("%s%d.%d.bin", prefix, i+1, j) //temporary because using hoons py file
				genofs[i][j] = gwas.NewGenoFileStream(filename, uint64(rows), uint64(blockSizes[j]), true)
			}
		}

		for i := 0; i < reg.K; i++ {
			for j := 0; j < reg.B; j++ {
				cols := foldSizes[i]
				filenameReg := fmt.Sprintf("%s%d.%d.bin", prefix, i+1, j)
				filename := fmt.Sprintf("%s_transpose%d.%d.bin", prefix, i+1, j)

				if _, err := os.Stat(filename); errors.Is(err, os.ErrNotExist) {
					log.LLvl1(time.Now(), "transposing")
					gwas.TransposeMatrixFile(filenameReg, cols, blockSizes[j], filename)
				}

				genoTfs[i][j] = gwas.NewGenoFileStream(filename, uint64(blockSizes[j]), uint64(cols), true)
			}
		}

		reg.general.SetGenoBlocks(genofs)
		reg.general.SetGenoTBlocks(genoTfs)
	}
}
func (reg *REGENIE) PrecomputeVals(step1 bool) {
	reg.LoadGFS(step1)
	network := reg.general.GetNetwork()
	cps := reg.general.GetCPS()
	Z_i := reg.general.GetCov()
	mpcObj := reg.general.GetMPC()
	pid := mpcObj.GetPid()
	var blockSizes, foldSizes []int
	if pid > 0 {
		blockSizes = reg.general.GetGenoBlockSizes()
		foldSizes = reg.general.GetGenoFoldSizes()
	}
	rtype := mpcObj.GetRType().Zero()

	fname := "mean_dosage.txt"

	var meanDosFloat []float64
	if reg.general.FileExistsForAll(mpcObj, reg.general.CacheFile(fname)) {
		if pid > 0 {
			meanDosDense := gwas.LoadMatrixFromFile(reg.general.CacheFile(fname), '\t')
			meanDosFloat = mat.DenseCopyOf(meanDosDense.T()).RawRowView(0)
		}
		log.LLvl1("Found mean dos")
	} else {
		mpcPar := reg.general.GetParallelMPC()
		numer := mpc_core.InitRVec(rtype, reg.M)
		denom := mpc_core.InitRVec(rtype, reg.M)

		if pid > 0 {
			numInds := uint64(reg.N)
			ac, _, miss := gwas.ReadGenoStatsFromFile(reg.general.GetConfig().GenoCountFile, reg.M)

			for i := range numer {
				numer[i] = rtype.FromUint64(uint64(ac[1][i]))
				denom[i] = rtype.FromUint64(numInds - uint64(miss[i]))
			}
		}

		meanDos := mpcPar.Divide(numer, denom)
		meanDos = mpcPar.RevealSymVec(meanDos)
		if pid > 0 {
			meanDosFloat = meanDos.ToFloat(mpcPar[0].GetFracBits())
			gwas.SaveFloatVectorToFile(reg.general.CacheFile(fname), meanDosFloat)
		}
	}

	if pid > 0 {
		genoBlockReplace := make([][]float64, reg.B)
		for i := 0; i < len(genoBlockReplace); i++ {
			genoBlockReplace[i] = make([]float64, blockSizes[i])
		}
		
		log.LLvl1("reg.B", reg.B)
		idx := 0
		for i := 0; i < reg.B; i++ {
			log.LLvl1("--i", i, len(genoBlockReplace[i]), len(meanDosFloat))
			for j := 0; j < len(genoBlockReplace[i]); j++ {
				log.LLvl1("----j", j, idx)
				genoBlockReplace[i][j] = meanDosFloat[idx]
				idx += 1
			}
		}

		// Create file streams for geno block files
		for i := 0; i < reg.K; i++ {
			for j := 0; j < reg.B; j++ {
				reg.general.GetGeno()[i][j].SetColMissingReplace(genoBlockReplace[j])
			}
		}
		for i := 0; i < reg.K; i++ {
			for j := 0; j < reg.B; j++ {
				reg.general.GetGenoT()[i][j].SetRowMissingReplace(genoBlockReplace[j])
			}
		}
	}

	var ZTZ crypto.CipherMatrix
	if pid > 0 {
		genoInvSig := reg.ComputeStdInv()
		invsig_k, _ := splitFolds(genoInvSig, blockSizes)
		invStdEncrBlock := make([]crypto.CipherVector, reg.B)
		for i := 0; i < len(invsig_k); i++ {
			invStdEncrBlock[i] = crypto.EncryptDense(cps, invsig_k[i])[0]
		}
		reg.invStdEncrBlock = invStdEncrBlock
		reg.invStd = mat.Col(nil, 0, genoInvSig)

		//compute Z_folds
		reg.Z_folds, reg.Z_kminus = splitFolds(Z_i, foldSizes)
		ZCache := make([]crypto.PlainMatrix, reg.K)
		ZkminusCache := make([]crypto.PlainMatrix, reg.K)

		for i := 0; i < reg.K; i++ {
			ZCache[i] = crypto.EncodeDense(cps, reg.Z_folds[i])
			ZkminusCache[i] = crypto.EncodeDense(cps, reg.Z_kminus[i])
		}
		reg.ZCache = ZCache
		reg.ZkminusCache = ZkminusCache
		//compute XTZ
		ScaledX_iTZ_iAll := make([]*mat.Dense, reg.B)

		nproc := runtime.GOMAXPROCS(0)

		jobChannels := make([]chan int, nproc)
		for i := range jobChannels {
			jobChannels[i] = make(chan int, 32)
		}

		// Dispatcher
		go func() {
			for b := 0; b < reg.B; b++ {
				jobChannels[b%nproc] <- b
			}

			for _, c := range jobChannels {
				close(c)
			}
		}()

		// Workers
		var workerGroup sync.WaitGroup
		for thread := 0; thread < nproc; thread++ {
			workerGroup.Add(1)
			go func(thread int) {
				defer workerGroup.Done()

				for b := range jobChannels[thread] {
					var sumXTZ *mat.Dense
					for k := 0; k < reg.K; k++ {
						geno_k_b := reg.general.GetGeno()[k][b]
						ScaledX_iTZ_i_k_b := gwas.MatTClear(geno_k_b, reg.Z_folds[k])
						ScaledX_iTZ_i_k_b.Scale(1/math.Sqrt(float64(reg.AllN)), ScaledX_iTZ_i_k_b)
						if k == 0 {
							sumXTZ = ScaledX_iTZ_i_k_b
						} else {
							sumXTZ.Add(sumXTZ, ScaledX_iTZ_i_k_b)
						}
					}
					ScaledX_iTZ_iAll[b] = sumXTZ
				}

			}(thread)
		}
		workerGroup.Wait()
		ScaledXTZcache := make([]crypto.PlainMatrix, reg.B)
		ScaledXTZ := make([]crypto.CipherMatrix, reg.B)
		ScaledXTZSS := make([]mpc_core.RMat, reg.B)

		// Dispatcher
		nproc = reg.general.GetMpcMainThreads()

		jobChannels = make([]chan int, nproc)
		for i := range jobChannels {
			jobChannels[i] = make(chan int, 32)
		}
		go func() {
			for b := 0; b < reg.B; b++ {
				jobChannels[b%nproc] <- b
			}

			for _, c := range jobChannels {
				close(c)
			}
		}()

		// Workers
		for thread := 0; thread < nproc; thread++ {
			workerGroup.Add(1)
			go func(thread int, cps *crypto.CryptoParams, mpcObj *mpc.MPC) {
				defer workerGroup.Done()
				for i := range jobChannels[thread] {
					slice := ScaledX_iTZ_iAll[i]
					ScaledXTZcache[i] = crypto.EncodeDense(cps, ScaledX_iTZ_iAll[i])
					sliceEncr := crypto.EncryptPlaintextMatrix(cps, ScaledXTZcache[i])
					ScaledXTZ[i] = mpcObj.Network.AggregateCMat(cps, sliceEncr)
					ScaledXTZSS[i] = make(mpc_core.RMat, reg.C)
					for j := range ScaledXTZSS[i] {
						vec := mat.VecDenseCopyOf(slice.ColView(j)).RawVector().Data
						ScaledXTZSS[i][j] = mpc_core.FloatToRVec(mpcObj.GetRType(), vec, mpcObj.GetFracBits())
					}
				}
			}(thread, cps.GetThread(thread), reg.general.GetParallelMPC()[thread])
		}
		workerGroup.Wait()

		reg.ScaledXTZ = ScaledXTZ
		reg.ScaledXTZSS = ScaledXTZSS
		reg.ScaledX_iTZ_i = ScaledX_iTZ_iAll
		reg.ScaledX_iTZ_iCache = ScaledXTZcache

		//compute ZTZinv
		Z_iTZ_i := mat.NewDense(reg.C, reg.C, nil)
		Z_iTZ_i.Mul(Z_i.T(), Z_i)

		encrZ_iTZ_i := crypto.EncryptDense(cps, Z_iTZ_i)
		ZTZ = network.AggregateCMat(cps, encrZ_iTZ_i)
		reg.ScaledZTZ = CMatMultConst(cps, ZTZ, 1/math.Sqrt(float64(reg.AllN)))
	}

	ss := mpc_core.InitRMat(rtype.Zero(), reg.C, reg.C)
	if pid > 0 {
		ss = mpcObj.CMatToSS(cps, mpcObj.GetRType(), ZTZ, -1, reg.C, 1, reg.C)
	}
	sqrtScale := rtype.FromFloat64(1/math.Sqrt(float64(reg.AllN)), mpcObj.GetFracBits())

	ss.MulScalar(sqrtScale)
	ss = mpcObj.TruncMat(ss, mpcObj.GetDataBits(), mpcObj.GetFracBits())

	ssInv, ssInvSqrt := mpcObj.MatrixInverseSymPos(ss)

	ssInv2 := ssInv.Copy()
	sqrtNScale := rtype.FromFloat64(math.Sqrt(float64(reg.AllN)), mpcObj.GetFracBits())
	ssInv2.MulScalar(sqrtNScale)
	ssInv2 = mpcObj.TruncMat(ssInv2, mpcObj.GetDataBits(), mpcObj.GetFracBits())

	// Scale up ssInvSqrt by sqrt(sqrt(N)) so that it is sqrt(N) times the real inverse sqrt
	doubleSqrtScale := rtype.FromFloat64(math.Sqrt(math.Sqrt(float64(reg.AllN))), mpcObj.GetFracBits())

	ssInvSqrt.MulScalar(doubleSqrtScale)
	ssInvSqrt = mpcObj.TruncMat(ssInvSqrt, mpcObj.GetDataBits(), mpcObj.GetFracBits())

	if pid > 0 {
		ScaledZTZinvTwice := mpcObj.SSToCMat(cps, ssInv2)
		reg.ScaledZTZinvTwice = ScaledZTZinvTwice
		reg.ScaledZTZinvSS = ssInv
		reg.ScaledZTZinvSqrtSS = ssInvSqrt
	}

	reg.ComputeTildePhenoFolds()
}
func (reg *REGENIE) ComputeTildePhenoFolds() {

	cps := reg.general.GetCPS()
	mpcObj := reg.general.GetMPC()
	rtype := mpcObj.GetRType()
	pid := mpcObj.GetPid()

	Z := reg.general.GetCov()
	y := reg.general.GetPheno()

	var tilde_pheno crypto.CipherVector
	tilde_pheno_k := make([]crypto.CipherVector, reg.K)
	tilde_pheno_kminus := make([]crypto.CipherVector, reg.K)

	var ZTySS mpc_core.RVec
	var R mpc_core.RMat
	if pid > 0 {
		ZTy := mat.NewDense(reg.C, 1, nil)
		ZTy.Mul(Z.T(), y)

		ZTySS = mpc_core.FloatToRVec(rtype, ZTy.RawMatrix().Data, mpcObj.GetFracBits())
		R = reg.ScaledZTZinvSS
	} else {
		ZTySS = mpc_core.InitRVec(rtype.Zero(), reg.C)
		R = mpc_core.InitRMat(rtype.Zero(), reg.C, reg.C)
	}
	ZZInvZTySS := mpcObj.SSMultMat(mpc_core.RMat{ZTySS}, R)
	ZZInvZTySS = mpcObj.TruncMat(ZZInvZTySS, mpcObj.GetDataBits(), mpcObj.GetFracBits())
	ZZInvZTy := mpcObj.SSToCVec(cps, ZZInvZTySS[0])

	var ySqSumSS mpc_core.RElem
	var Py crypto.CipherVector
	var ySqSum *ckks.Ciphertext
	Py_folds := make([]crypto.CipherVector, reg.K)
	Py_kminus := make([]crypto.CipherVector, reg.K)
	if pid > 0 {

		y_folds, y_kminus := splitFolds(y, reg.general.GetGenoFoldSizes())
		ZPlain := crypto.EncodeDense(cps, Z)
		Cy := gwas.CPMatMult5Parallel(cps, ZPlain, crypto.CipherMatrix{ZZInvZTy})[0]
		Cy = crypto.CMultConst(cps, Cy, -1/math.Sqrt(float64(reg.AllN)), true)
		yPlain, _ := crypto.EncodeFloatVector(cps, y.RawMatrix().Data)
		Py = crypto.CPAdd(cps, Cy, yPlain)

		for k := 0; k < reg.K; k++ {
			Cy_fold := gwas.CPMatMult5Parallel(cps, reg.ZCache[k], crypto.CipherMatrix{ZZInvZTy})[0]
			Cy_fold = crypto.CMultConst(cps, Cy_fold, -1/math.Sqrt(float64(reg.AllN)), true)
			yPlain_fold, _ := crypto.EncodeFloatVector(cps, y_folds[k].RawMatrix().Data)
			Py_folds[k] = crypto.CPAdd(cps, Cy_fold, yPlain_fold)

			ZPlain_kminus := crypto.EncodeDense(cps, reg.Z_kminus[k])
			Cy_kminus := gwas.CPMatMult5Parallel(cps, ZPlain_kminus, crypto.CipherMatrix{ZZInvZTy})[0]
			Cy_kminus = crypto.CMultConst(cps, Cy_kminus, -1/math.Sqrt(float64(reg.AllN)), true)
			yPlain_kminus, _ := crypto.EncodeFloatVector(cps, y_kminus[k].RawMatrix().Data)
			Py_kminus[k] = crypto.CPAdd(cps, Cy_kminus, yPlain_kminus)
		}
		ySqSum = crypto.InnerProd(cps, Py, Py)
		ySqSum = mpcObj.Network.AggregateCText(cps, ySqSum)
		ySqSum = crypto.CMultConst(cps, crypto.CipherVector{ySqSum}, 1/float64(reg.AllN), false)[0]
	} else {
		ySqSumSS = rtype.Zero()
	}

	ySqSumSS = mpcObj.CiphertextToSS(cps, rtype, ySqSum, -1, 1)[0]
	_, sigmaInvSS := mpcObj.SqrtAndSqrtInverse(mpc_core.RVec{ySqSumSS})
	sigmaInv := mpcObj.SStoCiphertext(cps, sigmaInvSS)
	sigmaInv = crypto.Rebalance(cps, sigmaInv)

	if pid > 0 {
		tilde_pheno = crypto.CMultScalar(cps, Py, sigmaInv)
		for k := 0; k < reg.K; k++ {
			tilde_pheno_k[k] = crypto.CMultScalar(cps, Py_folds[k], sigmaInv)
			tilde_pheno_kminus[k] = crypto.CMultScalar(cps, Py_kminus[k], sigmaInv)
		}
	}

	mpcObj.Network.BootstrapVecAll(cps, tilde_pheno)
	mpcObj.Network.BootstrapMatAll(cps, tilde_pheno_k)
	mpcObj.Network.BootstrapMatAll(cps, tilde_pheno_kminus)
	reg.tilde_pheno = tilde_pheno
	reg.tilde_pheno_k = tilde_pheno_k
}
