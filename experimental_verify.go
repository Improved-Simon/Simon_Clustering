package main

import (
	"crypto/rand"
	"encoding/binary"
	"flag"
	"fmt"
	"log"
	"math"
	"math/big"
	"os"
	"runtime"
	"sync"
	"time"
)

// ===================================================================================
// 1. INTERFACES AND GENERAL STRUCTURES
// ===================================================================================

type SimonCipher interface {
	BlockSize() int
	KeySize() int
	Encrypt(pt, k []byte) (y, x uint32)
	Rounds() int
}

type DifferentialCharacteristic struct {
	Name               string
	CipherVersion      string
	Epochs             int
	Rounds             int
	Delta              []uint32
	ExpectedOutputDiff []uint32
}

var characteristics = map[string]DifferentialCharacteristic{
	"64-10r-A": {
		Name: "SIMON64/128 10-round A", CipherVersion: "simon64", Epochs: 1 << 28, Rounds: 10,
		Delta: []uint32{0x40000004, 0x1}, ExpectedOutputDiff: []uint32{0x40000440, 0x1000},
	},
	"64-10r-B": {
		Name: "SIMON64/128 10-round B", CipherVersion: "simon64", Epochs: 1 << 27, Rounds: 10,
		Delta: []uint32{0x1880, 0x440}, ExpectedOutputDiff: []uint32{0x10100, 0x44040},
	},
	"64-10r-C": {
		Name: "SIMON64/128 10-round C", CipherVersion: "simon64", Epochs: 1 << 26, Rounds: 10,
		Delta: []uint32{0x222000, 0x80000}, ExpectedOutputDiff: []uint32{0x222000, 0x808000},
	},
	"64-10r-D": {
		Name: "SIMON64/128 10-round D", CipherVersion: "simon64", Epochs: 1 << 28, Rounds: 10,
		Delta: []uint32{0x11000000, 0x4000000}, ExpectedOutputDiff: []uint32{0x1000011, 0x40},
	},
	"48-13r-example": {
		Name: "SIMON48/96 8-round example", CipherVersion: "simon48", Epochs: 1 << 30, Rounds: 13,
		Delta: []uint32{0x000001, 0x000000}, ExpectedOutputDiff: []uint32{0x000001, 0x000000},
	},
	"32-13r-example": {
		Name: "SIMON32/64 13-round example", CipherVersion: "simon32", Epochs: 1 << 35, Rounds: 13,
		Delta: []uint32{0x0, 0x1}, ExpectedOutputDiff: []uint32{0x1, 0x0},
	},
}

var zSequences = map[string]uint64{
	"simon32": parseBitString("111110100010010101100001110011010000001100011011101011011001001"),
	"simon48": parseBitString("11011011100100101001111000001110101000010001100010110011111101"),
	"simon64": 0xfc2ce51207a635db,
}

func parseBitString(s string) uint64 { i, _ := new(big.Int).SetString(s, 2); return i.Uint64() }

// ===================================================================================
// 2. SIMON32/64 IMPLEMENTATION
// ===================================================================================

type Simon32Cipher struct{ rounds int }

func (c *Simon32Cipher) BlockSize() int { return 4 }
func (c *Simon32Cipher) KeySize() int   { return 8 }
func (c *Simon32Cipher) Rounds() int    { return c.rounds }
func (c *Simon32Cipher) Encrypt(pt, k []byte) (y, x uint32) {
	ptWords := make([]uint16, 2)
	kWords := make([]uint16, 4)
	bytesToWords16(pt, ptWords)
	bytesToWords16(k, kWords)
	rk := c.keySchedule(kWords)
	y16, x16 := c.encrypt(ptWords, rk)
	return uint32(y16), uint32(x16)
}
func (c *Simon32Cipher) keySchedule(K []uint16) []uint16 {
	rk := make([]uint16, c.rounds)
	const constC = uint16(0xfffc)
	z := zSequences["simon32"]
	copy(rk, K)
	for i := 4; i < c.rounds; i++ {
		tmp := rotr16(rk[i-1], 3) ^ rk[i-3]
		tmp ^= rotr16(rk[i-1], 4)
		tmp ^= rotr16(tmp, 1)
		rk[i] = constC ^ uint16(z&1) ^ rk[i-4] ^ tmp
		z >>= 1
	}
	return rk
}
func (c *Simon32Cipher) encrypt(Pt, rk []uint16) (y, x uint16) {
	y, x = Pt[0], Pt[1]
	numRounds := c.rounds

	// Execute rounds in pairs
	for i := 0; i < numRounds/2; i++ {
		roundIdx := i * 2
		y ^= f16(x) ^ rk[roundIdx]
		x ^= f16(y) ^ rk[roundIdx+1]
	}

	// If the number of rounds is odd, apply the final round
	if numRounds%2 != 0 {
		y ^= f16(x) ^ rk[numRounds-1]
	}

	return y, x
}
func rotl16(x uint16, r uint) uint16 { return (x << r) | (x >> (16 - r)) }
func rotr16(x uint16, r uint) uint16 { return (x >> r) | (x << (16 - r)) }
func f16(x uint16) uint16            { return (rotl16(x, 1) & rotl16(x, 8)) ^ rotl16(x, 2) }
func bytesToWords16(b []byte, w []uint16) {
	for i := range w {
		w[i] = binary.LittleEndian.Uint16(b[i*2:])
	}
}

// ===================================================================================
// 3. SIMON48/96 IMPLEMENTATION
// ===================================================================================

type Simon48Cipher struct{ rounds int }

func (c *Simon48Cipher) BlockSize() int { return 6 }
func (c *Simon48Cipher) KeySize() int   { return 12 }
func (c *Simon48Cipher) Rounds() int    { return c.rounds }
func (c *Simon48Cipher) Encrypt(pt, k []byte) (y, x uint32) {
	ptWords := make([]uint32, 2)
	kWords := make([]uint32, 4)
	bytesToWords24(pt, ptWords)
	bytesToWords24(k, kWords)
	rk := c.keySchedule(kWords)
	return c.encrypt(ptWords, rk)
}
func (c *Simon48Cipher) keySchedule(K []uint32) []uint32 {
	rk := make([]uint32, c.rounds)
	const constC = uint32(0xfffffc)
	z := zSequences["simon48"]
	copy(rk, K)
	for i := 4; i < c.rounds; i++ {
		tmp := rotr24(rk[i-1], 3) ^ rk[i-3]
		tmp ^= rotr24(rk[i-1], 4)
		tmp ^= rotr24(tmp, 1)
		rk[i] = constC ^ (uint32(z & 1)) ^ rk[i-4] ^ tmp
		z >>= 1
	}
	return rk
}
func (c *Simon48Cipher) encrypt(Pt, rk []uint32) (y, x uint32) {
	y, x = Pt[0], Pt[1]
	numRounds := c.rounds

	// Execute rounds in pairs
	for i := 0; i < numRounds/2; i++ {
		roundIdx := i * 2
		y ^= f24(x) ^ rk[roundIdx]
		x ^= f24(y) ^ rk[roundIdx+1]
	}

	// If the number of rounds is odd, apply the final round
	if numRounds%2 != 0 {
		y ^= f24(x) ^ rk[numRounds-1]
	}

	return y, x
}
func rotl24(x uint32, r uint) uint32 { return (((x << r) | (x >> (24 - r))) & 0xffffff) }
func rotr24(x uint32, r uint) uint32 { return (((x >> r) | (x << (24 - r))) & 0xffffff) }
func f24(x uint32) uint32            { return ((rotl24(x, 1) & rotl24(x, 8)) ^ rotl24(x, 2)) }
func bytesToWords24(b []byte, w []uint32) {
	for i := range w {
		w[i] = uint32(b[i*3+0]) | uint32(b[i*3+1])<<8 | uint32(b[i*3+2])<<16
	}
}

// ===================================================================================
// 4. SIMON64/128 IMPLEMENTATION
// ===================================================================================

type Simon64Cipher struct{ rounds int }

func (c *Simon64Cipher) BlockSize() int { return 8 }
func (c *Simon64Cipher) KeySize() int   { return 16 }
func (c *Simon64Cipher) Rounds() int    { return c.rounds }
func (c *Simon64Cipher) Encrypt(pt, k []byte) (y, x uint32) {
	ptWords := make([]uint32, 2)
	kWords := make([]uint32, 4)
	bytesToWords32(pt, ptWords)
	bytesToWords32(k, kWords)
	rk := c.keySchedule(kWords)
	return c.encrypt(ptWords, rk)
}
func (c *Simon64Cipher) keySchedule(K []uint32) []uint32 {
	rk := make([]uint32, c.rounds)
	const constC = uint32(0xfffffffc)
	z := zSequences["simon64"]
	copy(rk, K)
	for i := 4; i < c.rounds; i++ {
		tmp := rotr32(rk[i-1], 3) ^ rk[i-3]
		tmp ^= rotr32(rk[i-1], 4)
		tmp ^= rotr32(tmp, 1)
		rk[i] = constC ^ uint32(z&1) ^ rk[i-4] ^ tmp
		z >>= 1
	}
	return rk
}
func (c *Simon64Cipher) encrypt(Pt, rk []uint32) (y, x uint32) {
	y, x = Pt[0], Pt[1]
	numRounds := c.rounds

	// Execute rounds in pairs
	for i := 0; i < numRounds/2; i++ {
		roundIdx := i * 2
		y ^= f32(x) ^ rk[roundIdx]
		x ^= f32(y) ^ rk[roundIdx+1]
	}

	// If the number of rounds is odd, apply the final round
	if numRounds%2 != 0 {
		y ^= f32(x) ^ rk[numRounds-1]
	}

	return y, x
}
func rotl32(x uint32, r uint) uint32 { return (x << r) | (x >> (32 - r)) }
func rotr32(x uint32, r uint) uint32 { return (x >> r) | (x << (32 - r)) }
func f32(x uint32) uint32            { return (rotl32(x, 1) & rotl32(x, 8)) ^ rotl32(x, 2) }
func bytesToWords32(b []byte, w []uint32) {
	for i := range w {
		w[i] = binary.LittleEndian.Uint32(b[i*4:])
	}
}

// ===================================================================================
// 5. CONCURRENT EXPERIMENT FRAMEWORK
// ===================================================================================

func deltaToBytes(delta []uint32, cipher SimonCipher) []byte {
	b := make([]byte, cipher.BlockSize())
	switch cipher.(type) {
	case *Simon32Cipher:
		binary.LittleEndian.PutUint16(b[0:], uint16(delta[0]))
		binary.LittleEndian.PutUint16(b[2:], uint16(delta[1]))
	case *Simon48Cipher:
		w0, w1 := delta[0], delta[1]
		b[0], b[1], b[2] = byte(w0), byte(w0>>8), byte(w0>>16)
		b[3], b[4], b[5] = byte(w1), byte(w1>>8), byte(w1>>16)
	case *Simon64Cipher:
		binary.LittleEndian.PutUint32(b[0:], delta[0])
		binary.LittleEndian.PutUint32(b[4:], delta[1])
	default:
		log.Fatalf("deltaToBytes: unknown cipher type %T", cipher)
	}
	return b
}

func runExperimentWorker(
	epochsToRun int,
	cipher SimonCipher,
	delta []byte,
	ans []uint32,
	resultsChan chan<- int,
	wg *sync.WaitGroup,
) {
	defer wg.Done()
	localNum := 0
	pt1 := make([]byte, cipher.BlockSize())
	pt2 := make([]byte, cipher.BlockSize())
	k := make([]byte, cipher.KeySize())
	for i := 0; i < epochsToRun; i++ {
		if _, err := rand.Read(pt1); err != nil {
			log.Fatalf("worker: failed to generate random pt1: %v", err)
		}
		if _, err := rand.Read(k); err != nil {
			log.Fatalf("worker: failed to generate random k: %v", err)
		}
		for j := range pt1 {
			pt2[j] = pt1[j] ^ delta[j]
		}
		y1, x1 := cipher.Encrypt(pt1, k)
		y2, x2 := cipher.Encrypt(pt2, k)
		if ((y1 ^ y2) == ans[0]) && ((x1 ^ x2) == ans[1]) {
			localNum++
		}
	}
	resultsChan <- localNum
}

func main() {
	charName := flag.String("char", "64-10r-D", "The name of the differential characteristic to test")
	numWorkers := flag.Int("workers", runtime.NumCPU(), "Number of concurrent workers")
	listChars := flag.Bool("list", false, "List all available characteristics and exit")
	flag.Parse()

	if *listChars {
		fmt.Println("Available differential characteristics:")
		for name, char := range characteristics {
			fmt.Printf("  - %s:\n", name)
			fmt.Printf("    Cipher: %s, Rounds: %d, Epochs: 2^%.0f\n", char.CipherVersion, char.Rounds, math.Log2(float64(char.Epochs)))
			fmt.Printf("    Input Delta:  [0x%x, 0x%x]\n", char.Delta[0], char.Delta[1])
			fmt.Printf("    Output Diff:  [0x%x, 0x%x]\n", char.ExpectedOutputDiff[0], char.ExpectedOutputDiff[1])
		}
		os.Exit(0)
	}

	selectedChar, ok := characteristics[*charName]
	if !ok {
		log.Fatalf("Error: Characteristic '%s' not found. Use -list to see options.", *charName)
	}

	var cipher SimonCipher
	switch selectedChar.CipherVersion {
	case "simon32":
		cipher = &Simon32Cipher{rounds: selectedChar.Rounds}
	case "simon48":
		cipher = &Simon48Cipher{rounds: selectedChar.Rounds}
	case "simon64":
		cipher = &Simon64Cipher{rounds: selectedChar.Rounds}
	default:
		log.Fatalf("Error: Unknown cipher version '%s'", selectedChar.CipherVersion)
	}

	epochs := selectedChar.Epochs
	deltaBytes := deltaToBytes(selectedChar.Delta, cipher)
	ans := selectedChar.ExpectedOutputDiff

	fmt.Printf("Starting experiment: %s\n", selectedChar.Name)
	fmt.Printf("Using %d worker(s) on %d CPU(s)\n", *numWorkers, runtime.NumCPU())
	fmt.Printf("Total iterations (Epochs): 2^%.0f (%d)\n", math.Log2(float64(epochs)), epochs)
	fmt.Printf("Input Delta: [0x%x, 0x%x]\n", selectedChar.Delta[0], selectedChar.Delta[1])
	fmt.Printf("Expected Output Diff: [0x%x, 0x%x]\n", ans[0], ans[1])
	fmt.Println("-------------------------------------------------")

	startTime := time.Now()
	var wg sync.WaitGroup
	resultsChan := make(chan int, *numWorkers)
	epochsPerWorker := epochs / *numWorkers
	remainingEpochs := epochs % *numWorkers

	for i := 0; i < *numWorkers; i++ {
		wg.Add(1)
		epochsForThisWorker := epochsPerWorker
		if i == *numWorkers-1 {
			epochsForThisWorker += remainingEpochs
		}
		go runExperimentWorker(epochsForThisWorker, cipher, deltaBytes, ans, resultsChan, &wg)
	}

	go func() {
		wg.Wait()
		close(resultsChan)
	}()

	totalNum := 0
	for num := range resultsChan {
		totalNum += num
	}

	duration := time.Since(startTime)
	fmt.Printf("\nExperiment finished in: %v\n", duration)
	fmt.Printf("Total matches found (num): %d\n", totalNum)

	if totalNum > 0 {
		probability := float64(totalNum) / float64(epochs)
		log2Prob := math.Log2(probability)
		fmt.Printf("Empirical probability: %e\n", probability)
		fmt.Printf("Log2(Probability): %.4f\n", log2Prob)
	} else {
		fmt.Println("No matches found.")
	}
}
