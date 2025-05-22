package lmm

import (
	"fmt"
	"math"

	"go.dedis.ch/onet/v3/log"

	"strconv"
	"time"

	"github.com/hhcho/sfgwas-lmm/mpc"

	mpc_core "github.com/hhcho/mpc-core"
	"github.com/hhcho/sfgwas-lmm/gwas"
	"gonum.org/v1/gonum/mat"

	"github.com/hhcho/sfgwas-lmm/crypto"
	"github.com/ldsec/lattigo/v2/ckks"
)

//DUMMY FUNCTION
func (reg *REGENIE) RidgeRegressionCGD(cps *crypto.CryptoParams, mpcObj *mpc.MPC, lambdas []float64, fold, block int) crypto.CipherMatrix {
	var res crypto.CipherMatrix
	return res
}

func (reg *REGENIE) RidgeRegressionADMMWood(cps *crypto.CryptoParams, mpcObjs *mpc.ParallelMPC, lambdas []float64, Z_iTZ_iZTZinv crypto.CipherMatrix, ATA_k *mat.Dense, ScaledX_iTZ_i *mat.Dense, XTy crypto.CipherVector, fold, block int, calcGTicketsChannel chan int) crypto.CipherMatrix {
	start := time.Now()
	mpcObj := (*mpcObjs)[0]
	pid := mpcObj.GetPid()
	blockSize, runningSum := reg.GetBlockSizeInfo(block)
	fmt.Println("blockSize", blockSize)
	invStdSlice := make([]float64, blockSize)
	rho := 5000.0
	if pid > 0 {
		rho = float64(reg.AllN/reg.K) * 2
		copy(invStdSlice, reg.invStd[runningSum:runningSum+blockSize])
	}
	log.LLvl1(time.Now().Format(time.StampMilli), "Block", block, "Fold", fold, "Rho: ", rho, "ADMM regression compute local XTy,", time.Since(start))
	var ScaledXTZAll crypto.CipherMatrix
	var ScaledXTZCache crypto.PlainMatrix
	if pid > 0 {
		ScaledXTZAll = reg.ScaledXTZ[block]
		ScaledXTZCache = reg.ScaledX_iTZ_iCache[block]
	}

	params := &ADMMParams{
		general:        reg.general,
		cps:            cps,
		mpcObj:         mpcObj,
		mpcObjs:        mpcObjs,
		ATA:            ATA_k,
		ScaledX_iTZ_i:  ScaledX_iTZ_i,
		ScaledXTZAll:   ScaledXTZAll,
		ScaledXTZCache: ScaledXTZCache,
		Z_iTZ_iZTZinv:  Z_iTZ_iZTZinv,
		ScaledZTZ:      reg.ScaledZTZ,
		b:              XTy,
		invStd:         invStdSlice,
		rho:            rho,
		allN:           reg.AllN,
		C:              reg.C,
		T:              blockSize,
	}
	admm := NewADMM_Wood(params, block, fold, calcGTicketsChannel) // 26 min
	log.LLvl1(time.Now().Format(time.StampMilli), "Block", block, "Fold", fold, "ADMM regression make new ADMM: ", time.Since(start))
	predictors := make(crypto.CipherMatrix, len(lambdas))
	var lastRes crypto.CipherVector
	mpcObj.AssertSync()
	for i := 0; i < len(lambdas); i++ {
		predictors[i] = admm.Run(lambdas[i], math.Sqrt(lambdas[i])/rho, lastRes, reg.num_iteration_lvl0[i], block, fold) // change to scale
		predictors[i] = mpcObjs.GetNetworks().BootstrapVecAll(cps, predictors[i])
		if pid > 0 {
			for p := 1; p < mpcObj.GetNParty(); p++ {
				predictor_name := reg.generateBlockCache("predictor_"+strconv.Itoa(int(lambdas[i])), p, fold, block)
				gwas.SaveMatrixToFile(cps, mpcObj, crypto.CipherMatrix{predictors[i]}, blockSize, p, predictor_name)
			}
			lastRes = predictors[i]
		}
		mpcObj.AssertSync()
	}
	log.LLvl1(time.Now().Format(time.StampMilli), "Block", block, "Fold", fold, "RidgeRegressionADMMWood time: ", time.Since(start))
	return predictors
}

func (reg *REGENIE) RidgeRegressionLevel1(cps *crypto.CryptoParams, mpcObjs *mpc.ParallelMPC, Wkarr []crypto.CipherMatrix, mean, stdinv, WTy crypto.CipherVector, k int, lambdas []float64) crypto.CipherMatrix {
	start := time.Now()
	mpcObj := (*mpcObjs)[0]
	foldSizes := reg.general.GetGenoFoldSizes()
	pid := mpcObj.GetPid()
	out := make(crypto.CipherMatrix, len(lambdas))
	lastRes := CZeros(cps, reg.B*reg.R)

	for j := 0; j < len(lambdas); j++ {
		lazyMultFn := func(x crypto.CipherVector) crypto.CipherVector {
			sum := CZeros(cps, reg.B*reg.R)
			for i := 0; i < reg.K; i++ {
				if i == k {
					continue
				}
				Wx := reg.MultLazyCMatParallel(cps, mpcObjs, Wkarr[i], crypto.CipherMatrix{x}, mean, stdinv, foldSizes[k])[0]
				WxScaled := crypto.CMultConst(cps, Wx, 1/math.Sqrt(float64(reg.AllN)), false)
				WxScaled = mpcObjs.GetNetworks().BootstrapVecAll(cps, WxScaled)
				last := reg.MultLazyCMatTParallel(cps, mpcObjs, Wkarr[i], crypto.CipherMatrix{WxScaled}, mean, stdinv, foldSizes[k])[0]
				last = crypto.CMultConst(cps, last, 1/math.Sqrt(lambdas[j]), false) // to ensure that values don't blow up
				sum = crypto.CAdd(cps, sum, last)
			}
			newScale := math.Sqrt(lambdas[j]) / math.Sqrt(float64(reg.AllN)) // scale by sqrt(lambda[j])
			scaledX := crypto.CMultConst(cps, x, newScale, false)            //+ lambdaI
			sum = crypto.CAdd(cps, sum, scaledX)
			return sum
		}
		fmt.Println()
		log.LLvl1("Running for lambda ", lambdas[j], "index", j)
		fname := "etas" + strconv.Itoa(pid) + "_" + strconv.Itoa(j) + "_" + strconv.Itoa(k) + ".txt"
		if reg.general.FileExistsForAll(mpcObj, reg.general.OutFile(fname)) {
			log.LLvl1(time.Now().Format(time.StampMilli), "Output file found, skipping:", fname)
			if pid > 0 {
				out[j] = gwas.LoadCacheFromFile(cps, reg.general.OutFile(fname))[0]
			}
		} else {
			out[j], _, _ = ConjGradSolveCipherVec(cps, mpcObjs, WTy, reg.B*reg.R, lastRes, lazyMultFn, reg.num_iteration_lvl1[j], reg.lvl1_refresh_rate)
		}

		if pid > 0 {
			lastRes = out[j]
		}
		if pid > 0 {
			for p := 1; p < mpcObj.GetNParty(); p++ {
				gwas.SaveMatrixToFile(cps, mpcObj, crypto.CipherMatrix{out[j]}, reg.B*reg.R, p, reg.general.OutFile("etas"+strconv.Itoa(p)+"_"+strconv.Itoa(j)+"_"+strconv.Itoa(k)+".txt"))
				if p == pid {
					crypto.SaveCipherMatrixToFile(cps, crypto.CipherMatrix{out[j]}, reg.general.OutFile("cipher_eta"+strconv.Itoa(p)+"_"+strconv.Itoa(j)+"_"+strconv.Itoa(k)+".txt"))
				}
			}
		}
	}
	log.LLvl1(time.Now().Format(time.StampMilli), "RidgeRegressionLevel1 time: ", time.Now().Sub(start).String())

	return out
}

type LazyMult func(x crypto.CipherVector) crypto.CipherVector

//make into batches?
func ConjGradSolveCipherVec(cps *crypto.CryptoParams, mpcObjs *mpc.ParallelMPC, b crypto.CipherVector, b_len int, initialX crypto.CipherVector, LazyMult LazyMult, max_iter int, refresh_rate int) (crypto.CipherVector, int, bool) {
	useDummyBoot := false
	mpcObj := (*mpcObjs)[0]
	pid := mpcObjs.GetPid()
	//peeks at residual together after 3 iterations to see if there is progress
	RESPEEK := 10
	RESTHRES := 1e-3
	x := CZeros(cps, b_len)
	var p, r crypto.CipherVector
	if pid > 0 {
		if initialX == nil {
			r = crypto.CSub(cps, b, x)
			p = crypto.CopyEncryptedVector(r)
		} else {
			var initP crypto.CipherVector
			x = crypto.CopyEncryptedVector(initialX)
			initP = LazyMult(initialX)
			log.LLvl1(time.Now().Format(time.StampMilli), "init P")
			matPrint(cipherToNetworkMat(cps, mpcObj, crypto.CipherMatrix{initP}, b_len).T())
			if useDummyBoot {
				initP = mpcObj.Network.DummyBootstrapVecAll(cps, initP)
			} else {
				initP = mpcObjs.GetNetworks().BootstrapVecAll(cps, initP)
			}
			r = crypto.CSub(cps, b, initP)
			p = crypto.CopyEncryptedVector(r)
		}
	}

	var newX crypto.CipherVector
	k := 0
	for k = 0; k < max_iter; k++ {

		var Ap crypto.CipherVector
		if pid > 0 {
			if useDummyBoot {
				p = mpcObj.Network.DummyBootstrapVecAll(cps, p)
			} else {
				p = mpcObjs.GetNetworks().BootstrapVecAll(cps, p)
			}
			Ap = LazyMult(p)
		}
		Ap_scaled := crypto.CMultConst(cps, Ap, 1.0/math.Sqrt(float64(b_len)), false)
		numSS := CipherVectorInnerProdSS(cps, mpcObj, r, b_len)
		denomSS := CipherVectorInnerProd2SS(cps, mpcObj, p, Ap_scaled, b_len)

		var alpha *ckks.Ciphertext
		sqrtSS, sqrtInvSS := mpcObj.SqrtAndSqrtInverse(mpc_core.RVec{numSS, denomSS})
		fracSS := mpcObj.SSMultElemVec(mpc_core.RVec{sqrtSS[0]}, mpc_core.RVec{sqrtInvSS[1]})
		fracSS = mpcObj.TruncVec(fracSS, mpcObj.GetDataBits(), mpcObj.GetFracBits())

		fracSS = mpcObj.SSMultElemVec(fracSS, fracSS)
		fracSS = mpcObj.TruncVec(fracSS, mpcObj.GetDataBits(), mpcObj.GetFracBits())

		if pid > 0 {
			fracSS.MulScalar(mpcObj.GetRType().FromFloat64(1.0/math.Sqrt(float64(b_len)), mpcObj.GetFracBits()))
		}
		fracSS = mpcObj.TruncVec(fracSS, mpcObj.GetDataBits(), mpcObj.GetFracBits())

		alpha = mpcObj.SStoCiphertext(cps, fracSS)
		alpha = crypto.Rebalance(cps, alpha)

		var newR crypto.CipherVector

		if pid > 0 {
			scaledP := crypto.CMultScalar(cps, p, alpha)
			newX = crypto.CAdd(cps, x, scaledP)
			if useDummyBoot {
				newX = mpcObj.Network.DummyBootstrapVecAll(cps, newX)
			} else {
				newX = mpcObjs.GetNetworks().BootstrapVecAll(cps, newX)
			}
			if k != 0 && k%refresh_rate == 0 {
				Ax := LazyMult(newX)
				newR = crypto.CSub(cps, b, Ax)
			} else {
				scaledAp := crypto.CMultScalar(cps, Ap, alpha)
				newR = crypto.CSub(cps, r, scaledAp)
			}

			if useDummyBoot {
				newR = mpcObj.Network.DummyBootstrapVecAll(cps, newR)
			} else {
				newR = mpcObjs.GetNetworks().BootstrapVecAll(cps, newR)
			}
		}

		newRDotSS := CipherVectorInnerProdSS(cps, mpcObj, newR, b_len)
		oldRDotSS := CipherVectorInnerProdSS(cps, mpcObj, r, b_len)

		if pid > 0 {
			r = newR
			x = newX
		}

		// Check if time to peek
		signal := CONTINUE

		if k > 0 {
			log.LLvl1("Initiating error peek, k", k, "RESPEEK", RESPEEK)

			//decrypt and break if residual is low enough
			var newRDotPlain float64
			if pid > 0 {
				newRDotPlain = mpcObj.RevealSym(newRDotSS).Float64(mpcObj.GetFracBits())
				log.LLvl1("ERROR PEEK: ", newRDotPlain)

				pv := mpcObj.Network.CollectiveDecryptVec(cps, newR, 1)
				fv := crypto.DecodeFloatVector(cps, pv)
				S := 0.0
				for ii := range fv {
					S += fv[ii] * fv[ii]
				}
				log.LLvl1("Expected: ", S)

				if newRDotPlain < -10 {
					log.LLvl1("DEBUG: NEGATIVE ERROR!")
					log.Fatal()
				}
			}

			if pid == mpcObj.GetHubPid() {
				if math.Abs(newRDotPlain) <= RESTHRES {
					signal = STOP
				} else {
					signal = CONTINUE
				}

				for other := 0; other < mpcObj.GetNParty(); other++ {
					if other != pid {
						mpcObj.Network.SendInt(signal, other)
					}
				}
			} else { // All other parties including PID = 0
				signal = mpcObj.Network.ReceiveInt(mpcObj.GetHubPid())
			}
		}

		if signal == STOP {
			log.LLvl1("Converged! Terminating early at iter", k)
			return x, k, true
		}

		var beta *ckks.Ciphertext
		sqrtSSBeta, sqrtInvSSBeta := mpcObj.SqrtAndSqrtInverse(mpc_core.RVec{newRDotSS, oldRDotSS})

		fracSSBeta := mpcObj.SSMultElemVec(mpc_core.RVec{sqrtSSBeta[0]}, mpc_core.RVec{sqrtInvSSBeta[1]})
		fracSSBeta = mpcObj.TruncVec(fracSSBeta, mpcObj.GetDataBits(), mpcObj.GetFracBits())

		fracSSBeta = mpcObj.SSMultElemVec(fracSSBeta, fracSSBeta)
		fracSSBeta = mpcObj.TruncVec(fracSSBeta, mpcObj.GetDataBits(), mpcObj.GetFracBits())

		beta = mpcObj.SStoCiphertext(cps, fracSSBeta)
		beta = crypto.Rebalance(cps, beta)

		if pid > 0 {
			scaledBeta := crypto.CMultScalar(cps, p, beta)
			p = crypto.CAdd(cps, newR, scaledBeta)
		}
	}
	return x, k, false
}
