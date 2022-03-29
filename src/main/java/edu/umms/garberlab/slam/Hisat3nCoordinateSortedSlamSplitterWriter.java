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

import java.io.File;
import java.io.Serializable;
import java.util.HashMap;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMProgramRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;

public class Hisat3nCoordinateSortedSlamSplitterWriter extends Hisat3nSplitterWriter{
	private SAMFileWriter unconvertedWriter;
	private SAMFileWriter convertedWriter;
	
	private SamReader pairQueryReader;
	
	
	private String unconvertedFilePath;
	private String unconvertedFileIdxPath;
	private String convertedFilePath;
	private String convertedFileIdxPath;
	
	//private SAMFileHeader header;
	private FragmentCache cache;


	public Hisat3nCoordinateSortedSlamSplitterWriter(File outDirFile, String outPrefix, SamReader reader, boolean asBAM)  {
		this.unconvertedFilePath = outDirFile.getAbsolutePath() + "/" + outPrefix + "_unconverted.bam";
		this.unconvertedFileIdxPath = outDirFile.getAbsolutePath() + "/" + outPrefix + "_unconverted.bai";
		this.convertedFilePath = outDirFile.getAbsolutePath() + "/" + outPrefix + "_converted.bam";
		this.convertedFileIdxPath = outDirFile.getAbsolutePath()+ "/" + outPrefix + "_converted.bai";
		
		SAMFileHeader header = reader.getFileHeader();
		this.pairQueryReader = reader;
		
		SAMProgramRecord pg = new SAMProgramRecord(SlamSplitter.PROGRAM_NAME);
		header.addProgramRecord(pg);
		
		SAMFileWriterFactory factory = new SAMFileWriterFactory();
		if(asBAM) {
			unconvertedWriter = factory.makeBAMWriter(header, false, new File(unconvertedFilePath));
			convertedWriter = factory.makeBAMWriter(header, false, new File(convertedFilePath));
                    
		} else {
			unconvertedWriter = factory.makeSAMWriter(header, false, new File(unconvertedFilePath));
			convertedWriter = factory.makeSAMWriter(header, false, new File(convertedFilePath));
		}
		
		cache = new FragmentCache();
			
	}
	
	public Hisat3nCoordinateSortedSlamSplitterWriter(File outDirFile, String outPrefix, SAMFileHeader header, boolean asBAM)  {
		this.unconvertedFilePath = outDirFile.getAbsolutePath() + "/" + outPrefix + "_unconverted.bam";
		this.unconvertedFileIdxPath = outDirFile.getAbsolutePath() + "/" + outPrefix + "_unconverted.bai";
		this.convertedFilePath = outDirFile.getAbsolutePath() + "/" + outPrefix + "_converted.bam";
		this.convertedFileIdxPath = outDirFile.getAbsolutePath()+ "/" + outPrefix + "_converted.bai";
		
		
		SAMProgramRecord pg = new SAMProgramRecord(SlamSplitter.PROGRAM_NAME);
		header.addProgramRecord(pg);
		
		SAMFileWriterFactory factory = new SAMFileWriterFactory();
		if(asBAM) {
			unconvertedWriter = factory.makeBAMWriter(header, false, new File(unconvertedFilePath));
			convertedWriter = factory.makeBAMWriter(header, false, new File(convertedFilePath));
                    
		} else {
			unconvertedWriter = factory.makeSAMWriter(header, false, new File(unconvertedFilePath));
			convertedWriter = factory.makeSAMWriter(header, false, new File(convertedFilePath));
		}
		
			
	}
	
	
	public void close() {
		unconvertedWriter.close();
		convertedWriter.close();
		
	}
	
	public void write(SAMRecord samRecord) {
		Hisat3nAlignedFragment fragment = null;
		if(isPairedEnd()) {
			fragment = cache.getFragment(samRecord.getReadName());
		}
		if(fragment == null) {
			fragment = new Hisat3nAlignedFragment(new Hisat3nSAMRecord(samRecord), pairQueryReader);
		}
		
		write(fragment);

	}

	public void write(Hisat3nAlignedFragment hisatAlignmentFragment) {
		
		
		if(!isPairedEnd() ||  hisatAlignmentFragment.isUnpaired()) {
			write(hisatAlignmentFragment.getPair1());
		} else {
			Hisat3nSAMRecord upstreamRecord = hisatAlignmentFragment.getUpstreamRecord();
			Hisat3nSAMRecord downstreamRecord = hisatAlignmentFragment.getDownstreamRecord();
			if(cache.hasFragment(downstreamRecord.getSAMRecord().getReadName())) {
				write(downstreamRecord, hisatAlignmentFragment.getConvertedBases()>0);
				cache.removeFragment(downstreamRecord.getSAMRecord().getReadName());
			} else {
				cache.addFragment(hisatAlignmentFragment);
				write(upstreamRecord, hisatAlignmentFragment.getConvertedBases()>0);
			}
		}
		
	}
	
	private void write(Hisat3nSAMRecord record) { 
		write (record, record.countConvertedBases()>0);
	}
	
	private void write(Hisat3nSAMRecord record, boolean hasConvertedBases) {
		if(hasConvertedBases) {
			if(revertConversion()) {
				record.revertConvertedBases();
			}
			convertedWriter.addAlignment(record.getSAMRecord());
		} else {
			unconvertedWriter.addAlignment(record.getSAMRecord());
		}
	}
	

}

class FragmentCache implements Serializable {

	private static final long serialVersionUID = 2806066906990741041L;

	// Change to HashTable if/when making it a threaded application
	HashMap<String, Hisat3nAlignedFragment> writeCache = new HashMap<String, Hisat3nAlignedFragment>(); 
	
	void addFragment(Hisat3nAlignedFragment fragment) {
		if(!fragment.isUnpaired()) {
			writeCache.put(fragment.getPair1().getSAMRecord().getReadName(), fragment);
		}
	}

	public Hisat3nAlignedFragment getFragment(String readName) {
		return writeCache.get(readName);
	}

	public boolean hasFragment(String readName) {
		return writeCache.containsKey(readName);
	}
	
	public void removeFragment(String readName) {
		writeCache.remove(readName);
	}
	
}