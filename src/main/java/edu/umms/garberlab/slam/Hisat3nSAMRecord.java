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

import java.io.Serializable;
import java.util.Iterator;
import java.util.List;
import java.util.regex.MatchResult;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;

public class Hisat3nSAMRecord implements Serializable{
	/**
	 * 
	 */
	private static final long serialVersionUID = 7599727388327500299L;
	//private static final Pattern forwardConversionMDPattern = Pattern.compile("T[0-9]+|[0-9]+T");
	//private static final Pattern reverseConversionMDPattern = Pattern.compile("A[0-9]+|[0-9]+A");
	private static final Pattern MDBlockPattern = Pattern.compile("[0-9]+[A|C|G|T]|[0-9]+\\^[[A|C|G|T]]+");
	
	public final static int MIN_INTRON_LENGTH = 20;
	public final static int MIN_MATCH_LENTH = 15;
	
	private String mdTag;
	private int nhTag;
	private int nmTag;
	private int asTag;
	private char yzTag;
	private int yfTag;
	
	private SAMRecord samRecord;
	
	private Iterator<CigarElement> cigarElntIterator; //This is highly thread unsafe. But could not find a different way
	private int cigarIdx;

	public Hisat3nSAMRecord(SAMRecord samRecord) {
		this.samRecord = samRecord;
		if(!samRecord.getReadUnmappedFlag()) {
	        mdTag = (String) samRecord.getAttribute("MD");
	        nhTag = (int) samRecord.getAttribute("NH");
	        nmTag = (int) samRecord.getAttribute("NM");
	        asTag = (int) samRecord.getAttribute("AS");
	        yzTag = (char) samRecord.getAttribute("YZ");
	        yfTag = (int) samRecord.getAttribute("Yf");
		}
		cigarElntIterator = samRecord.getCigar().iterator();
		cigarIdx = 0;
	}


	public String getMdTag() {
		return mdTag;
	}


	public void setMdTag(String mdTag) {
		this.mdTag = mdTag;
	}


	public int getNhTag() {
		return nhTag;
	}


	public void setNhTag(int nhTag) {
		this.nhTag = nhTag;
	}


	public int getNmTag() {
		return nmTag;
	}


	public void setNmTag(int nmTag) {
		this.nmTag = nmTag;
	}


	public int getAsTag() {
		return asTag;
	}


	public void setAsTag(int asTag) {
		this.asTag = asTag;
	}


	public char getYzTag() {
		return yzTag;
	}


	public void setYzTag(char yzTag) {
		this.yzTag = yzTag;
	}


	public int getYfTag() {
		return yfTag;
	}


	public void setYfTag(int yfTag) {
		this.yfTag = yfTag;
	}


	public boolean isSpliced() {
        Cigar cigar = samRecord.getCigar();	
        List<CigarElement> elements = cigar.getCigarElements();
        
        boolean hasAMatch = false;
        boolean hasGap    = false;
        
        for (CigarElement element : elements) {
     	   if(CigarOperator.N.equals(element.getOperator()) && element.getLength() > MIN_INTRON_LENGTH) {
     		   hasGap = true;
     	   }
     	   if(CigarOperator.M.equals(element.getOperator()) && element.getLength() > MIN_MATCH_LENTH) {
     		   hasAMatch = true;
     	   }
        }
        
		return hasAMatch && hasGap;
	}


	/**
	 * In the Hisat3n version that was used to implement this function
	 * the aligner reports T->C conversions on reads aligned to the (-) strand
	 * as converted. This is not right. Here we check and ensure that at least one of the
	 * converted bases is a T->C for + aligned reads or a G->A for - aligned reads
	 * 
	 * ACTUALLY AFTER RUNNING THIS VERSION I NOTICED THAT THE OBSERVATION ABOVE WAS NOT TRUE
	 * HISAT2-3N is not counting the wrong conversions! 
	 * 
	 * COMMENTING THIS OUT
	 *
	public int countConvertedBases() {
		int convertedBases = 0;
		if (getYfTag()>0) {
			char strand = getYzTag();
			String matchInfo = getMdTag();
			Matcher matcher = '+' == strand ? forwardConversionMDPattern.matcher(matchInfo)  : reverseConversionMDPattern.matcher(matchInfo); 
			while(matcher.find()) {
				convertedBases ++;
			}
			
		}
		
		return convertedBases;
	}
	*/
	
	public int countConvertedBases() {
		return getYfTag();
	}
	
	public float getFractionOfConvertedBases() {
		int convertedBases = countConvertedBases();
		return convertedBases > 0 ? convertedBases/getYzTag() : 0;
	}


	public void revertConvertedBases() {
		Cigar cigar = samRecord.getCigar();
		byte [] readBases = samRecord.getReadBases();
		String read =samRecord.getReadString();
		char strand = getYzTag();
		String matchInfo = getMdTag();
		char convertToBaseType = '+' == strand ? 'T' : 'A';
		char convertedToBase     = '+' == strand ? 'C' : 'G';
		byte convertedToBaseByte = (byte) convertedToBase;
		
		int idx = 0;
		this.cigarIdx = 0;
		//Looking for soft clipped bases
		CigarElement first = this.cigarElntIterator.next();
		if (first.getOperator().equals(CigarOperator.S)) {
			idx = idx + first.getLength();
			this.cigarIdx = idx + this.cigarElntIterator.next().getLength();	//I am assuming that only matches can follow clipped bases
		}
		
		
		// Now walk through the MD
		Matcher matcher = MDBlockPattern.matcher(matchInfo);
		short revertedBases = 0;
		while(matcher.find()) {
			String match = matcher.group();
			if(match.contains("^")) { 
				String [] deletionInfo = match.split("\\^");
				int numOfPreviousMatches = Integer.parseInt(deletionInfo[0]);
				//int lengthOfInsertion = deletionInfo[1].length();
				idx = incrementReadIdx(idx, numOfPreviousMatches); //Ignore the reference insertion, it does no affect the read idx
			} else {
				int blockLength = Integer.parseInt(match.substring(0, match.length() - 1));
				idx = incrementReadIdx(idx , blockLength);
				char mismatchCharacter = match.charAt(match.length() - 1);
				if(mismatchCharacter == convertToBaseType  && readBases[idx] == convertedToBaseByte) {
					readBases[idx] = (byte) mismatchCharacter;
					revertedBases++;
				}
				idx = incrementReadIdx(idx,1); //This one is for the actual base that was converted regardless of whether or not the base needed to be reverted
			}			
		}
		if(revertedBases != getYfTag()) {
			System.err.println(samRecord.format());
			throw new RuntimeException("BUG: Reverted "+revertedBases +" when there were "+countConvertedBases()+" bases converted");
		}
		samRecord.setReadBases(readBases);		
		String newread = samRecord.getReadString();
		
		
	}


	private int incrementReadIdx(int idx, 
			int numOfPreviousMatches) {
		
		idx = idx+numOfPreviousMatches;
		while(this.cigarIdx < idx && this.cigarElntIterator.hasNext()) {
			CigarElement elmnt = this.cigarElntIterator.next();
			CigarOperator elmntOp = elmnt.getOperator();

			//Need to handle hard clipped, X and other operators
			switch (elmntOp) {
			case I:
				idx = idx + elmnt.getLength();
			case M:
				this.cigarIdx = cigarIdx + elmnt.getLength();
				break;
			case N:
			case D:
			default:
				break;
			}			
		}
		
		return idx;
	}


	public SAMRecord getSAMRecord() {
		return this.samRecord;
	}
	
	public boolean isMapped() {
		return ! this.samRecord.getReadUnmappedFlag();
	}
	
	
	
	
}