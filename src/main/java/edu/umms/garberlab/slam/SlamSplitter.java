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
import java.util.List;
import java.util.Map.Entry;
import java.util.Set;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMProgramRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

/**
 * This program seeks to sort converted and uncoverted reads from an alignment file 
 * resulting from hisat-3n by writing out reads.
 * 
 * Reads can be optionally unconverted to help downstream processing
 * 
 * To identify and manipulate converted reads, the program relies on HISAT-3N specific attributes
 * as described in http://daehwankimlab.github.io/hisat2/hisat-3n/:
 * 
 * Yf:i:<N>: Number of conversions are detected in the read.
 * Zf:i:<N>: Number of un-converted bases are detected in the read. Yf + Zf = total number of bases which can be converted in the read sequence.
 * YZ:A:<A>: The value + or – indicate the read is mapped to REF-3N (+) or REF-RC-3N (-).
 * 
 * @author Manuel Garber
 *
 */

public class SlamSplitter {
	
	public static String PROGRAM_NAME = "SlamSlpitter";
	
	
	public static void main(String[] args) {
		final Options options = new Options();
		final Option help = new Option("h","help", false, "print this message");
		final Option writeSam = new Option("s","sam", false, "Write SAM (BAM is default)");
		
		final Option inputBamOpt = Option.builder("inputBam")
				.argName("bam")
				.hasArg()
				.required(true)
				.desc("Required - Input Bam File")
				.build();
		
		final Option outputDirectoryOpt = Option.builder("outDir")
				.argName("outDir")
				.hasArg()
				.required(false)
				.desc("Output directory if other than current")
				.build();

		final Option outputPrefixOpt = Option.builder("outPrefix")
				.argName("outPrefix")
				.hasArg()
				.required(false)
				.desc("Required - prefix for the splitted output files")
				.build();
		
		final Option gentotype = Option.builder("genotypeVCF")
				.argName("genotypeVCF")
				.hasArg()
				.desc("VCF file containing genotype information to handle T/C variants")
				.build();
		
		options.addOption(inputBamOpt)
		.addOption(outputDirectoryOpt)
		.addOption(outputPrefixOpt)
		.addOption(gentotype)
		.addOption(help)
		.addOption(writeSam); 
		
		CommandLine cmd;
		CommandLineParser clip = new DefaultParser();
		HelpFormatter helper = new HelpFormatter();
		
		long totalSplicedAlignments = 0;
		long totalAlignments        = 0;
		long totalConvertedReads    = 0;
		long totalUnmapped          = 0;
		try {
			cmd = clip.parse(options, args);
			if(cmd.hasOption(help)) {
				helper.printHelp("", options);
				System.exit(0);
			}
			
			
			String inBam = cmd.getOptionValue(inputBamOpt);
			File inBamFile = new File(inBam);
			if (!inBamFile.exists()) {
				System.err.println("The bam file path provided "+inBam+" does not exist");
				System.exit(1);
			}
			
			File outDirFile = new File(cmd.getOptionValue(outputDirectoryOpt, "."));
			if (!outDirFile.exists()) {
				System.err.println("The output directory  provided "+outDirFile.getAbsolutePath()+" does not exist");
				System.exit(1);
			}
			
			String [] inBamPathComponents = inBam.split("/");
			String outPrefix = cmd.getOptionValue(outputPrefixOpt, 
					inBamPathComponents[inBamPathComponents.length - 1].replaceFirst("\\.[^\\.]+$", ""));
			
			
			SamReader reader = SamReaderFactory.makeDefault().open(new File(inBam));
			
			SAMFileHeader header = reader.getFileHeader();
			
			Hisat3nSlamSplitterWriter sw = new Hisat3nSlamSplitterWriter(outDirFile, outPrefix, header, !cmd.hasOption(writeSam));
			sw.turnOnConversionReversion();
			
			Set<Entry<String, String>> attributes = header.getAttributes();
			
			for (Entry<String,String> att : attributes) {
				System.out.println("Atribute key " + att.getKey() + " and val "+ att.getValue());
			}
			
			List<SAMProgramRecord>programRecords = header.getProgramRecords();
			boolean isHisat3N = false;
			for (SAMProgramRecord pgRcrd : programRecords) {
				System.out.println("Program record id " + pgRcrd.getId() + " and cmd: "+ pgRcrd.getCommandLine());
				isHisat3N = isHisat3N || (pgRcrd.getId().contains("hisat2") && pgRcrd.getProgramVersion().contains("-3n"));  // Kind of an unfortunate hack
				
			}
			
			if (!isHisat3N) { System.err.println("Warning - this program is meant to get Hisat-3n output. It does not seem that this alignment was generated by it");}
			
			
	        for (final SAMRecord samRecord : reader) {
	        	totalAlignments++;
	        	Hisat3nSAMRecord hisatRecord = new Hisat3nSAMRecord(samRecord);
	        	
	        	if(hisatRecord.isMapped()) {
		        	if (hisatRecord.isSpliced()) {
		        		totalSplicedAlignments++;
		        	}
		        	
		        	if(hisatRecord.countConvertedBases() > 0) {
		        		totalConvertedReads++;
		        	}
	        	} else {
	        		totalUnmapped++;
	        	}

	        	sw.write(hisatRecord);
	           
	           //samRecord.
	           
	        }
	        
	        sw.close();
			
			
		} catch (ParseException e) {
            //System.out.println(e.getMessage());
            helper.printHelp("", options);
            System.exit(0);
        } 
		
		System.out.println("Total alignments: " + totalAlignments + 
				"\n\t\tslpliced: " + totalSplicedAlignments + 
				"\n\t\ttotal reads with conversion: " + totalConvertedReads +
				"\n\t\ttotal unmapped reads: " + totalUnmapped);
		
		
	}
	

}
