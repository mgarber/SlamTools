/*
 * The MIT License
 *
 * Copyright (c) 2009 The Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */

package edu.umms.garberlab.slam;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;

/**
 * @author mgarber
 * Represents the alginmnet of a fragment. Has pointers to both the first and second paris when such exist
 */
public class Hisat3nAlignedFragment {
	
	static final int DEFAULT_FRAGMENT_SIZE = 300;
	static final int DEFAULT_FRAGMENT_SIZE_SD = 300;
	
	private Hisat3nSAMRecord pair1;
	private Hisat3nSAMRecord pair2;
	
	private int convertedBases;
	private double libraryFragmentLengthMean;
	private double libraryFragmentLengthStdDev;
	
	public Hisat3nAlignedFragment(Hisat3nSAMRecord hisatRecord, SamReader pairQueryReader) {
		libraryFragmentLengthMean = DEFAULT_FRAGMENT_SIZE;
		libraryFragmentLengthStdDev = DEFAULT_FRAGMENT_SIZE_SD;
	 	if(hisatRecord.getSAMRecord().getMateAlignmentStart() > 0) {
	 		SAMRecord mateSamRecord = pairQueryReader.queryMate(hisatRecord.getSAMRecord());
	 		if(mateSamRecord == null) {
	 			System.err.println("ERROR: mate had a non-zero alignment start but queryMate returned null\n" + hisatRecord.getSAMRecord().format());
	 		}
	 		Hisat3nSAMRecord mate = new Hisat3nSAMRecord(mateSamRecord);
	 		if(mate.getSAMRecord().getFirstOfPairFlag()) {
	 			pair1 = mate;
	 			pair2 = hisatRecord;
	 		} else {
	 			pair1 = hisatRecord;
	 			pair2 = mate;
	 		}
	 		convertedBases = pair1.countConvertedBases()+pair2.countConvertedBases();
	 	} else {
	 		pair1 = hisatRecord;
	 		convertedBases = pair1.countConvertedBases();
	 	}
	}
	
	public Hisat3nAlignedFragment(SAMRecord rec1, SAMRecord rec2) {
		libraryFragmentLengthMean = DEFAULT_FRAGMENT_SIZE;
		libraryFragmentLengthStdDev = DEFAULT_FRAGMENT_SIZE_SD;
		Hisat3nSAMRecord hisatRec1 = new Hisat3nSAMRecord(rec1);
		Hisat3nSAMRecord hisatRec2 = new Hisat3nSAMRecord(rec2);

		if(rec1.getFirstOfPairFlag()) {
			pair1 = hisatRec1;
			pair2 = hisatRec2;
		} else {
			pair1 = hisatRec2;
			pair2 = hisatRec1;
		}
		
		convertedBases = pair1.countConvertedBases()+pair2.countConvertedBases();

	}

	public Hisat3nSAMRecord getPair1() {
		return pair1;
	}

	public Hisat3nSAMRecord getPair2() {
		return pair2;
	}
	
	public double getLibraryFragmentLengthMean() {
		return libraryFragmentLengthMean;
	}

	public void setLibraryFragmentLengthMean(double libraryFragmentLengthMean) {
		this.libraryFragmentLengthMean = libraryFragmentLengthMean;
	}

	public double getLibraryFragmentLengthStdDev() {
		return libraryFragmentLengthStdDev;
	}

	public void setLibraryFragmentLengthStdDev(double libraryFragmentLengthStdDev) {
		this.libraryFragmentLengthStdDev = libraryFragmentLengthStdDev;
	}

	public boolean isUnpaired() {
		return pair2 == null;
	}

	public Hisat3nSAMRecord getUpstreamRecord() {
		Hisat3nSAMRecord retVal = null;
		if(!isUnpaired()) {
			retVal =  pair1.getSAMRecord().getAlignmentStart() < pair2.getSAMRecord().getAlignmentStart() ? pair1 : pair2;
		} 
		return retVal;
	}
	
	public Hisat3nSAMRecord getDownstreamRecord() {
		Hisat3nSAMRecord retVal = null;
		if(!isUnpaired()) {
			retVal =  pair1.getSAMRecord().getAlignmentStart() < pair2.getSAMRecord().getAlignmentStart() ? pair2 : pair1;
		} 
		return retVal;
	}

	public int getConvertedBases() {
		return convertedBases;
	}
	
	public boolean isSpliced () {
		return pair1.isSpliced() || 
				(pair2 != null && pair2.isSpliced()) ||
				pair2 != null && 
					(getDownstreamRecord().getSAMRecord().getAlignmentEnd() - 
							getUpstreamRecord().getSAMRecord().getAlignmentStart() 
							> getLibraryFragmentLengthMean() + 2*getLibraryFragmentLengthStdDev());
							
	}
	
	public boolean isMapped() {
		return pair1.isMapped() ||
				(pair2 != null && pair2.isMapped());
	}
	
}
