/*
 * The MIT License
 *
 * Copyright (c) 2022 University of Massachusetts Medical School
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
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;

public class Hisat3nSAMRecord implements Serializable{
	/**
	 * 
	 */
	private static final long serialVersionUID = 7599727388327500299L;
	//private static final Pattern forwardConversionMDPattern = Pattern.compile("T[0-9]+|[0-9]+T");
	//private static final Pattern reverseConversionMDPattern = Pattern.compile("A[0-9]+|[0-9]+A");
	private static final Pattern MDBlockPattern = Pattern.compile("[0-9]+\\^[[A|C|G|T|N]]+|[0-9]+[A|C|G|T|N]?");
	
	public final static int MIN_INTRON_LENGTH = 20;
	public final static int MIN_MATCH_LENTH = 15;
	
	private String mdTag;
	private int nhTag;
	private int nmTag;
	private int asTag;
	private char yzTag;
	private int yfTag;
	
	private SAMRecord samRecord;
	
	//private Iterator<CigarElement> cigarElntIterator; //This is highly thread unsafe. But could not find a different way
	private int cigarElementArrayIdx;
	private int cigarIdx;
	
	private Map<Integer, Character> sequenceToReferenceMap ;

	public Hisat3nSAMRecord(SAMRecord samRecord) {
		this.samRecord = samRecord;
		sequenceToReferenceMap = new HashMap<Integer, Character> ();
		if(!samRecord.getReadUnmappedFlag()) {
	        mdTag = (String) samRecord.getAttribute("MD");
	        nhTag = (int) samRecord.getAttribute("NH");
	        nmTag = (int) samRecord.getAttribute("NM");
	        asTag = (int) samRecord.getAttribute("AS");
	        yzTag = (char) samRecord.getAttribute("YZ"); 
	        yfTag = (int) samRecord.getAttribute("Yf");
	        sequenceToReferenceMap = makeSequenceToReferenceMap(samRecord.getCigar(), mdTag);
		}
		//cigarElntIterator = samRecord.getCigar().iterator();
		cigarElementArrayIdx = 0;
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
		char strand = getYzTag();
		byte [] readBases = samRecord.getReadBases();
		char convertToBaseType = '+' == strand ? 'T' : 'A';
		char convertedToBase     = '+' == strand ? 'C' : 'G'; 
		byte convertedToBaseByte = (byte) convertedToBase;
		short revertedBases = 0;
		
		for(int pos : sequenceToReferenceMap.keySet()) {
			char mismatchCharacter = sequenceToReferenceMap.get(pos);
			if(mismatchCharacter == convertToBaseType  && readBases[pos] == convertedToBaseByte) {
				readBases[pos] = (byte) mismatchCharacter;
				revertedBases++;
			}
		}
		
		if(revertedBases != getYfTag()) {			
			System.err.println("BUG: - reverted bases mismatch\n"+samRecord.format());
			throw new RuntimeException("BUG: Reverted "+revertedBases +" when there were "+countConvertedBases()+" bases converted ");
		}
	}


	public void revertConvertedBases2() {
		Cigar cigar = samRecord.getCigar();
		if(samRecord.getReadUnmappedFlag()) {
			return;
		}
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
		//CigarElement first = this.cigarElntIterator.next();
		CigarElement first = this.samRecord.getCigar().getCigarElement(this.cigarElementArrayIdx);
		this.cigarElementArrayIdx = this.cigarElementArrayIdx + 1;
		if (first.getOperator().equals(CigarOperator.S)) {
			idx = idx + first.getLength();
			//I am assuming that only matches can follow clipped bases. The cigar idx only advances on matches not insertions.
			//this.cigarIdx = first.getLength() + this.cigarElntIterator.next().getLength();	
			this.cigarIdx = first.getLength() +  
					this.samRecord.getCigar().getCigarElement(this.cigarElementArrayIdx).getLength();
		} else {
			this.cigarIdx = first.getLength();
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
			
			System.err.println("BUG: - reverted bases mismatch\n"+samRecord.format());
			throw new RuntimeException("BUG: Reverted "+revertedBases +" when there were "+countConvertedBases()+" bases converted ");
		}
		samRecord.setReadBases(readBases);				
	}


	
	private int incrementReadIdx(int idx, 
			int numOfPreviousMatches) {
		
		idx = idx+numOfPreviousMatches;
		
		//this is tricky. if an I occurs we need to move the cigar index again. so it would move twice, hence the while
		//while(this.cigarIdx <= idx && this.cigarElntIterator.hasNext()) { 
		while(this.cigarIdx <= idx && cigarElementArrayIdx < getSAMRecord().getCigar().getCigarElements().size()) {
			//CigarElement elmnt = this.cigarElntIterator.next(); 
			CigarElement elmnt = samRecord.getCigar().getCigarElement(cigarElementArrayIdx);
			cigarElementArrayIdx=cigarElementArrayIdx+1;
		
			CigarOperator elmntOp = elmnt.getOperator();

			//Need to handle hard clipped, X and other operators
			switch (elmntOp) {
			case I:
				idx = idx + elmnt.getLength();
				break;
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
	
	

	protected static Map<Integer, Character> makeSequenceToReferenceMap (Cigar cigar, String mdTag) {
		Map<Integer, Character> queryMismatch = new HashMap<Integer, Character>();
		Iterator<CigarElement>  cigarElmntIt = cigar.iterator();
		Matcher matcher = MDBlockPattern.matcher(mdTag);
		//String [] reconstruction = new String[cigar.getReadLength(cigar.getCigarElements())];
		int queryIdx = 0;
		int cigarCummulativeWalk = 0;
		//I am assuming that only matches can follow clipped bases. The cigar idx only advances on matches not insertions.
		//this.cigarIdx = first.getLength() + this.cigarElntIterator.next().getLength();	
		CigarElement cigarElement = cigarElmntIt.next();

		if(cigarElement.getOperator().equals(CigarOperator.S) ) {
			//We advance all index references
			queryIdx =  cigarElement.getLength();
			cigarCummulativeWalk = cigarElement.getLength();
			for (int i = 0; i< queryIdx ; i++) {
				queryMismatch.put(i, '-');
			}
			cigarElement = cigarElmntIt.next();
		}

		//NOTE: because we dealt with soft clipped bases at the begginign. We *must* be in a match element.
		//I am assuming that the first cigar element (other than a soft clip) must be a match.
		cigarCummulativeWalk += cigarElement.getLength();
		CigarOperator co = cigarElement.getOperator(); 
		if(!co.equals(CigarOperator.M)) {
			throw new RuntimeException("Expected a M (match) element as the first cigar element after other than soft clip but got " + co.toString());
		}
		String mdElement = null;
		char mismatch = '*'; // use Z for no mismatch;
		while(matcher.find()) {
			mdElement = matcher.group();
			
			mismatch = '*';
			int blockLength = 0;
			if(matcher.hitEnd()) {
				blockLength = Integer.parseInt(mdElement);
			} else if (mdElement.contains("^")){
				String [] refInsertionInfo = mdElement.split("\\^");
				
				try {
					blockLength = Integer.parseInt(refInsertionInfo[0]);	
				} catch (java.lang.NumberFormatException e) {
					System.err.println("Error parsing md tag.\nMD: "+ mdTag +"\ncigar: " + cigar.toString());
					e.printStackTrace();
					System.exit(1)
;				}
				//mismatch =  mdElement.charAt(mdElement.length() - 1);
			} else {
				mismatch =  mdElement.charAt(mdElement.length() - 1);
				blockLength = Integer.parseInt(mdElement.substring(0, mdElement.length() - 1));
			}

			for (int i = 0; i < blockLength; i++) {
				if(queryIdx == cigarCummulativeWalk && cigarElmntIt.hasNext()) {
					cigarElement =  cigarElmntIt.next();
					co = cigarElement.getOperator();				


					if(co.equals(CigarOperator.D) || co.equals(CigarOperator.N)) {
						cigarElement = cigarElmntIt.next(); // After a reference insertion there must be another element;
						co = cigarElement.getOperator();
					}
					if(co.equals(CigarOperator.I)) {
						for (int j = 0; j < cigarElement.getLength(); j++) {
							queryMismatch.put(queryIdx, '-');
							queryIdx += 1;
							cigarCummulativeWalk += 1;
						}
						cigarElement = cigarElmntIt.next(); //After a query insertion there must be another element;
						co = cigarElement.getOperator();
					}
					cigarCummulativeWalk += cigarElement.getLength();
				}
				queryIdx += 1;
			}
			
			if(queryIdx == cigarCummulativeWalk && cigarElmntIt.hasNext()) {
				cigarElement =  cigarElmntIt.next();
				co = cigarElement.getOperator();				

				if(co.equals(CigarOperator.D) || co.equals(CigarOperator.N)) {
					cigarElement = cigarElmntIt.next(); // After a reference insertion there must be another element;
					co = cigarElement.getOperator();
				}

				if(co.equals(CigarOperator.I)) {
					for (int j = 0; j < cigarElement.getLength(); j++) {
						queryMismatch.put(queryIdx, '-');
						queryIdx += 1;
						cigarCummulativeWalk += 1;
					}
					cigarElement = cigarElmntIt.next(); //After a query insertion there must be another element;
					co = cigarElement.getOperator();
				}


				cigarCummulativeWalk += cigarElement.getLength();
			}
			if(! ('^' == mismatch || '*' == mismatch) ) {
				queryMismatch.put(queryIdx, mismatch);
				queryIdx += 1;
			}
			// Done updating Cigar position, now lets look at the MD tag


		}
		return queryMismatch;
	}
	
}