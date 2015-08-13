#! /nethome/sravishankar/local/bin/python3
#BISSNP pipeline version1
import os
import re
import sys
import glob
import config
import shutil
import logging
import argparse
import subprocess
from multiprocessing import Process

def trimGalore(fastq, output_dir, fastq_type, adaptor, assay, log_file,basename):
    """Execute TrimGalore, initial QC, determine quality statistics, perform quality based sequence trimming\n
    Input to module: Path to fastq file, output directory path, assay type(RRBS/WGBS), fastq type(single/paired), log file path"""
    logging.info('Running trimgalore')
    logger = open(log_file,'a')
    trim_dir = output_dir 
    fastqc_dir = output_dir 
    a1 = adaptor
    trim_cmd = [config.trim_galore,'-o',trim_dir,'--phred33', '--fastqc_args','\"--outdir='+fastqc_dir+'\"','-a',a1]
    if fastq_type == 'paired':
        fastq1 = fastq[0]
        fastq2 = fastq[1]
        a2 = adaptor
        trim_cmd += ['-a2',a2,'--paired',fastq1,fastq2]
    elif fastq_type == 'single':
        fastq1 = fastq[0]
        trim_cmd.append(fastq1)
    if assay == 'RRBS':
        trim_cmd += ['--rrbs']
    run_trimgalore = subprocess.Popen(trim_cmd,shell=False,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    trimgalore_screen = run_trimgalore.communicate()
    logger.write('Trimgalore STDOUT:\n'+trimgalore_screen[0].decode("utf-8")+'\nTrimgalore STDERR:\n'+trimgalore_screen[1].decode("utf-8")+'\n')
    logger.close()
    return(glob.glob(trim_dir+'/*f*q'))

def bwaMeth(trim_fastq, fastq_type, output_dir, threads, reference, log_file,basename):
    """Execute bwaMeth, bisulfite seq specific alignment program,\n
    Input to module: Path to trimmed fastq files, fastq file type, number of threads, reference file path, output directory path"""
    logging.info('Running bwaMeth')
    logger = open(log_file,'a')
    bwa_out = output_dir 
    bwa_cmd = ['/nethome/sravishankar/local/bin/python3',config.bwa_meth,'--prefix',bwa_out+'/'+basename,'--thread',threads,'--reference',reference]
    if fastq_type == 'single':
        bwa_cmd.append(trim_fastq[0])
    elif fastq_type == 'paired':
        bwa_cmd += trim_fastq
    print(' '.join(bwa_cmd))
    run_bwameth = subprocess.Popen(bwa_cmd,shell=False, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    bwameth_screen = run_bwameth.communicate()
    logger.write('BWAMeth STDOUT:\n'+bwameth_screen[0].decode("utf-8")+'\nBWAMeth STDERR:\n'+bwameth_screen[1].decode("utf-8")+'\n')
    logger.close()
    return(glob.glob(bwa_out+'/'+basename+'*bam')[0])

def addReadGroup(bwa_bam, output_dir, reference, log_file, tempdir,basename,mem):
    logging.info('Running AddOrReplaceReadGroup')
    logger = open(log_file,'a')
    arg_out = output_dir
    out_file = arg_out+'/'+basename+'_readGroup.bam'
    baseref = os.path.splitext(os.path.basename(reference))[0]
    arg_cmd = [config.java,'-Djava.io.tmpdir='+tempdir,mem,'-jar',config.arg,'INPUT='+bwa_bam,'OUTPUT='+out_file,'RGID='+basename,'RGLB='+basename,'RGPL=\"Illumina\"',
            'RGPU=\"HiSeq-2500\"','RGSM='+basename,'RGDS='+baseref,'RGCN=\"NYGenome\"','CREATE_INDEX=true']
    run_arg = subprocess.Popen(arg_cmd,shell=False,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    arg_screen = run_arg.communicate()
    logger.write('AddReadGroup STDOUT:\n'+arg_screen[0].decode("utf-8")+'\nAddReadGroup STDERR:\n'+arg_screen[1].decode("utf-8")+'\n')
    logger.close()
    return(out_file)

def mergeBam(bam_files, output_dir, log_file, basename, tempdir,mem):
    logging.info('Running MergeSamFiles')
    logger = open(log_file,'a')
    mbam_out = output_dir+'/bissnp_results'
    if not os.path.exists(mbam_out):
        os.mkdir(mbam_out)
    out_file = mbam_out + '/'+basename+'_finalbam.bam'
    comment = ','.join(bam_files)
    input_string = ['INPUT='+filename for filename in bam_files]
    mbam_cmd = [config.java, '-Djava.io.tmpdir='+tempdir,mem,'-jar',config.mbam,input_string,'OUTPUT='+out_file,'USE_THREADING=\'true\'','COMMENT='+comment]
    run_mbam = subprocess.Popen(mbam_cmd,shell=False,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    mbam_screen = run_mbam.communicate()
    logger.write('MergeBamFile STDOUT:\n'+mbam_screen[0].decode("utf-8")+'\nMergeBamFile STDERR:\n'+mbam_screen[1].decode("utf-8")+'\n')
    logger.close()
    return(out_file)

def bisRTC(final_bam, output_dir, reference, known_files, interval_file, threads, log_file,basename,mem):
    logging.info('Running BisSNP RealignerTarget')
    logger = open(log_file,'a')
    brtc_out = output_dir 
    known_string = ' '.join(['--known '+ filename for filename in known_files]).split(' ')
    if interval_file != None:
        brtc_cmd = [config.java,mem,'-jar',config.bissnp,'-T','BisulfiteRealignerTargetCreator','-R',reference, '-I',final_bam,
                '-L',interval_file,'-o',brtc_out+'/'+basename+'_indel_target_interval.intervals','-nt',threads] + known_string
    else:
        brtc_cmd = [config.java,mem,'-jar',config.bissnp,'-T','BisulfiteRealignerTargetCreator','-R',reference, '-I',final_bam,
                    '-o',brtc_out+'/'+basename+'_indel_target_interval.intervals','-nt',threads] + known_string
    print(' '.join(brtc_cmd))
    run_brtc = subprocess.Popen(brtc_cmd, shell=False, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    #run_brtc.wait()
    brtc_screen = run_brtc.communicate()
    logger.write('BisulfiteRealignerTargetCreator STDOUT:\n'+brtc_screen[0].decode("utf-8")+'\nBisulfiteRealignerTargetCreator STDERR:\n'+brtc_screen[1].decode("utf-8")+'\n')
    logger.close()
    return(brtc_out+'/'+basename+'_indel_target_interval.intervals')

def bisInRealign(final_bam, reference, brtc_out, known_files, output_dir, log_file,basename,mem):
    logging.info('Running BisSNP IndelRealigner')
    logger = open(log_file,'a')
    bir_out = output_dir 
    known_string = ' '.join(['-known ' + filename for filename in known_files]).split(' ') #temp fix for known vcf
    if not os.path.exists(bir_out):
        os.mkdir(bir_out)
    bir_cmd = [config.java,mem,'-jar', config.bissnp,'-T','BisulfiteIndelRealigner','-R',reference, '-I',final_bam,
            '-targetIntervals',brtc_out,'-cigar','-o',bir_out+'/'+basename+'_realigned.bam'] + known_string
    print(' '.join(bir_cmd))
    bir_run = subprocess.Popen(bir_cmd,shell=False,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    bir_screen = bir_run.communicate()
    logger.write('BisulfiteIndelRealigner STDOUT:\n'+bir_screen[0].decode("utf-8")+'\nBisulfiteIndelRealigner STDERR:\n'+bir_screen[1].decode("utf-8")+'\n')
    logger.close()
    return(bir_out+'/'+basename+'_realigned.bam')

def markDuplicate(bir_out, output_dir, log_file,basename,mem):
    logging.info('Running MarkDuplicate')
    logger = open(log_file,'a')
    md_out = output_dir
    outfile = md_out + '/' + basename +'_mdup.bam'
    metrics = md_out + '/' + basename +'_metric.txt'
    md_cmd = [config.java, mem,'-jar',config.md,'I='+bir_out,'O='+outfile,'METRICS_FILE='+metrics,'CREATE_INDEX=true','VALIDATION_STRINGENCY=SILENT']
    md_run = subprocess.Popen(md_cmd, shell=False, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    md_screen = md_run.communicate()
    logger.write('MarkDuplicate STDOUT:\n'+md_screen[0].decode("utf-8")+'\nMarkDuplicate STDERR:\n'+md_screen[1].decode("utf-8")+'\n')
    logger.close()
    return(outfile)

def countCovariant(dup_out, reference, dbsnp, threads, log_file,output_dir,suffix, basename, mem):
    logging.info('Running count covariant')
    logger = open(log_file,'a')
    cc_out = output_dir
    outfile = cc_out + '/' + basename +'_'+suffix+'_recal.csv'
    known_string = ' '.join(['-knownSites ' + filename for filename in dbsnp]).split(' ') #temp fix for known vcf
    cc_cmd = [config.java, mem,'-jar',config.bissnp,'-T','BisulfiteCountCovariates','-R',reference, '-I',dup_out,
        '-cov','ReadGroupCovariate','-cov','QualityScoreCovariate','-cov','CycleCovariate',
            '-recalFile',outfile,'-nt',threads] + known_string
    print (cc_cmd)
    cc_run = subprocess.Popen(cc_cmd, shell=False, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    cc_screen = cc_run.communicate()
    ana_dir = cc_out + '/'+basename+'_analyze_covariates_'+suffix
    os.mkdir(ana_dir)
    ana_cmd = [config.java, mem,'-jar','/nethome/sravishankar/tools/BisulfiteAnalyzeCovariates-0.69.jar','-recalFile',outfile,'-outputDir',ana_dir,'--ignoreQ','5','--max_quality_score','40']
    ana_run = subprocess.Popen(ana_cmd,shell=False,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    logger.write('BisulfiteCountCovariates STDOUT:\n'+cc_screen[0].decode("utf-8")+'\nBisulfiteCountCovariates STDERR:\n'+cc_screen[1].decode("utf-8")+'\n')
    logger.close()
    return(outfile)

def recalBQS(cc_out,dup_out,reference,output_dir,log_file,basename,threads,mem):
    logging.info('Running BQSR')
    logger = open(log_file,'a')
    bqsr_out = output_dir 
    outfile = bqsr_out + '/' + basename +'_recal.bam'
    bqsr_cmd = [config.java,mem,'-jar',config.bissnp,'-T','BisulfiteTableRecalibration','-R',reference, '-I',dup_out,'-o',outfile,'-recalFile',cc_out,'-maxQ','40','-nt',threads]
    bqsr_run = subprocess.Popen(bqsr_cmd,shell=False,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    bqsr_screen = bqsr_run.communicate()
    logger.write('BisulfiteTableRecalibration STDOUT:\n'+bqsr_screen[0].decode("utf-8")+'\nBisulfiteTableRecalibration STDERR:\n'+bqsr_screen[1].decode("utf-8")+'\n')
    logger.close()
    return(outfile)

def bisSNP(reference,bqsr_out,dbsnp,intervals,log_file,output_dir,basename,threads, mem):
    logging.info('Running bisSNP')
    logger = open(log_file,'a')
    bsnp_out = output_dir
    outfile_cpg = bsnp_out + '/' + basename + '_cpg.vcf'
    outfile_snp = bsnp_out + '/' + basename + '_snp.vcf'
    bsnp_cmd = [config.java,mem,'-jar',config.bissnp,'-R',reference,'-T','BisulfiteGenotyper',
            '-I',bqsr_out,'-D',dbsnp,'-vfn1',outfile_cpg,'-vfn2',outfile_snp,
            '-stand_call_conf','20','-stand_emit_conf','10','-mmq','20','-mbq','20','-nt',threads]  #'-L',intervals
    bsnp_run = subprocess.Popen(bsnp_cmd, shell=False, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    bsnp_screen = bsnp_run.communicate()
    logger.write('BisulfiteGenotyper STDOUT:\n'+bsnp_screen[0].decode("utf-8")+'\nBisulfiteGenotyper STDERR:\n'+bsnp_screen[1].decode("utf-8")+'\n')
    logger.close()
    return(outfile_cpg,outfile_snp)

if __name__ == '__main__':
    logging.basicConfig(format='%(asctime)s %(message)s',level=logging.INFO)
    logging.info('BisSNP pipeline started')
    parser = argparse.ArgumentParser(description='BISSNP pipeline')
    parser.add_argument('-f','--fastq',nargs='+',type=str,help='Fastq file path')
    parser.add_argument('-o','--output',type=str,dest='output_dir',help='Output directory path')
    parser.add_argument('-t','--type',type=str,choices=['single','paired'],help='Fastq file type')
    parser.add_argument('-a','--adaptor',type=str,default='AGATCGGAAGAGC',help='Adaptor sequence')
    parser.add_argument('-s','--assay',type=str,help='Assay type')
    parser.add_argument('-n','--threads',type=str,help='Number of threads',default='8')
    parser.add_argument('-r','--reference',type=str,help='Reference file path')
    parser.add_argument('-k','--known',type=str,nargs='+',help='Known indel files')
    parser.add_argument('-l','--intervals',type=str,help='Intervals file')
    parser.add_argument('-d','--dbsnp',type=str,help='dbSNP vcf file path')
    parser.add_argument('-m','--javamem',type=str,help='Java heap memory',default='10g')
    args = parser.parse_args()
    file_name = os.path.basename(args.fastq[0])
    basename = re.split('[_.]R1[.]+',file_name)[0]
    log_file = args.output_dir +'/'+basename+'_logfile.txt'
    xmx = '-Xmx'+args.javamem
    #Creating directory structure
    temp_dir = args.output_dir + '/temp_dir'
    qc_dir = args.output_dir + '/qc'
    analysis_dir = args.output_dir + '/analysis'
    if not os.path.exists(temp_dir):
        os.mkdir(temp_dir)
    if not os.path.exists(qc_dir):
        os.mkdir(qc_dir)
    if not os.path.exists(analysis_dir):
        os.mkdir(analysis_dir)
    trim_fastq = trimGalore(args.fastq,qc_dir,args.type,args.adaptor,args.assay,log_file,basename)
    bwa_bam = bwaMeth(trim_fastq,args.type, analysis_dir,args.threads,args.reference,log_file,basename)
    arg_bam = addReadGroup(bwa_bam, analysis_dir, args.reference, log_file, temp_dir,basename,xmx)
    #bam_files = glob.glob(bamfile_dir+'/*_readGroup.bam')
    #if len(bam_files) > 1:
    #    final_bam = mergeBam(bam_files, args.output_dir, log_file, basename, temp_dir)
    #elif len(bam_files) == 1:
    #    final_bam_dir = args.output_dir + '/bissnp_results'
    #    if not os.path.exists(final_bam_dir):
    #        os.mkdir(final_bam_dir)
    #    final_bam = final_bam_dir + '/' + basename +'_finalbam.bam'
    #    shutil.move(bam_files[0],final_bam)
    #    sam_cmd = ['samtools','index',final_bam]
    #    run_sam = subprocess.Popen(sam_cmd,shell=False)
    #    run_sam.wait()
    brtc_out = bisRTC(arg_bam, analysis_dir, args.reference, args.known, args.intervals, args.threads, log_file, basename,xmx)
    bir_out = bisInRealign(arg_bam, args.reference, brtc_out, args.known, analysis_dir, log_file,basename,xmx)
    dup_out = markDuplicate(bir_out,analysis_dir,log_file,basename,xmx)
    cc_out = countCovariant(dup_out, args.reference, args.known, args.threads,log_file,analysis_dir,'pre',basename,xmx)
    bqsr_out = recalBQS(cc_out,dup_out,args.reference,analysis_dir,log_file,basename,args.threads,xmx)
    cc_out = Process(target=countCovariant,args=(bqsr_out, args.reference, args.known, args.threads,log_file,analysis_dir,'post',basename,xmx))
    cc_out.start()
    bsnp_cpg,bsnp_snp = bisSNP(args.reference, bqsr_out, args.dbsnp, args.intervals, log_file,analysis_dir,basename,args.threads,xmx)
    cc_out.join()
    os.remove(arg_bam)
    os.remove(bir_out)
    os.remove(dup_out)
