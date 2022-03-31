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

import java.util.regex.Matcher;
import java.util.regex.Pattern;

import htsjdk.samtools.SAMRecord;

/**
 * @author mgarber
 *
 */
 abstract  class Hisat3nSplitterWriter {
	
	private boolean revertConversion;
	private boolean isPairedEnd;
	private boolean writeUnmapped;

	abstract public void close();
	abstract public void write(SAMRecord samRecord);
	abstract public void write(Hisat3nAlignedFragment hisatAlignmentFragment);
	
	
	private static final Pattern ReadNameWithSpace = Pattern.compile(" +.*");
		
	public void turnOnConversionReversion() {
		this.revertConversion = true;
	}
	
	public void turnOfConversionReversion() {
		this.revertConversion = false;
	}
	
	protected boolean revertConversion() {
		return this.revertConversion;
	}
	
	public void setPairedEnd(boolean pairedEnd) {
		this.isPairedEnd = pairedEnd;
	}

	public boolean isPairedEnd() {
		return this.isPairedEnd ;
	}
	public  void setWriteUnmapped(boolean writeUnmapped) {
		this.writeUnmapped = writeUnmapped;
	}
	public boolean isWriteUnmapped() {
		return writeUnmapped;
	}
	
	protected static boolean readNamesMatch(SAMRecord samRecord, SAMRecord samRecord2) {
		return readNamesMatch(samRecord.getReadName(), samRecord2.getReadName());
	}
	
	protected static boolean readNamesMatch(String samRecordReadName1, String samRecordReadName2) {
		
		if (samRecordReadName1.contains(" ") && samRecordReadName2.contains(" ")) {
			Matcher r1m = ReadNameWithSpace.matcher(samRecordReadName1);
			samRecordReadName1 = r1m.replaceAll("");
			
			Matcher r2m = ReadNameWithSpace.matcher(samRecordReadName2);
			samRecordReadName2 = r2m.replaceAll("");
			
		}
		return samRecordReadName1.equals(samRecordReadName2);
	}

}
