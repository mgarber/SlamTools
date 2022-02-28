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

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;

import htsjdk.samtools.BAMStreamWriter;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMFileWriterImpl;
import htsjdk.samtools.SAMProgramRecord;

public class Hisat3nSlamSplitterWriter {
	private SAMFileWriter unconvertedWriter;
	private SAMFileWriter convertedWriter;
	
	
	private String unconvertedFilePath;
	private String unconvertedFileIdxPath;
	private String convertedFilePath;
	private String convertedFileIdxPath;
	
	private SAMFileHeader header;
	
	private boolean revertConversion;


	public Hisat3nSlamSplitterWriter(File outDirFile, String outPrefix, SAMFileHeader header, boolean asBAM) {
		this.unconvertedFilePath = outDirFile.getAbsolutePath() + "/" + outPrefix + "_unconverted.bam";
		this.unconvertedFileIdxPath = outDirFile.getAbsolutePath() + "/" + outPrefix + "_unconverted.bai";
		this.convertedFilePath = outDirFile.getAbsolutePath() + "/" + outPrefix + "_converted.bam";
		this.convertedFileIdxPath = outDirFile.getAbsolutePath()+ "/" + outPrefix + "_converted.bai";
		
		this.header = header;
		
		SAMProgramRecord pg = new SAMProgramRecord(SlamSplitter.PROGRAM_NAME);
		this.header.addProgramRecord(pg);
		
		SAMFileWriterFactory factory = new SAMFileWriterFactory();
		if(asBAM) {
			unconvertedWriter = factory.makeBAMWriter(header, false, new File(unconvertedFilePath));
			convertedWriter = factory.makeBAMWriter(header, false, new File(convertedFilePath));
                    
		} else {
			unconvertedWriter = factory.makeSAMWriter(header, false, new File(unconvertedFilePath));
			convertedWriter = factory.makeSAMWriter(header, false, new File(convertedFilePath));
		}
			
	}
	
	public void turnOnConversionReversion() {
		this.revertConversion = true;
	}
	
	public void turnOfConversionReversion() {
		this.revertConversion = false;
	}
	

	public void close() {
		unconvertedWriter.close();
		convertedWriter.close();
		
	}

	public void write(Hisat3nSAMRecord hisatRecord) {
		if(hisatRecord.countConvertedBases()>0) {
			if(revertConversion) {
				hisatRecord.revertConvertedBases();
			}
			convertedWriter.addAlignment(hisatRecord.getSAMRecord());
		} else {
			unconvertedWriter.addAlignment(hisatRecord.getSAMRecord());
		}
		
	}
	

}
