#!/usr/bin/env python
import argparse
import sys, os
import pickle as p

class common:
    @staticmethod
    def system(cmd):
        print >> sys.stderr, cmd
        if os.system(cmd):
            print >> sys.stderr, "Error with command:", cmd
            sys.exit(-1)

    @staticmethod
    def create_path(path):
        if not os.path.exists(path): os.makedirs(path)
        return path

    @staticmethod
    def gbk2fasta(gbk, fasta):
        out = open(fasta, 'w')
        found_origin = False
        for line in open(gbk, 'r'):
            if line.startswith("LOCUS"):
                out.write('>' + line.split()[1] + '\n')
            elif line.startswith("ORIGIN"):
                found_origin = True
            elif found_origin:
                tokens = line.split()
                seq = "".join(tokens[1:]).upper()
                out.write(seq + '\n')
        common.assert_file(fasta)

    @staticmethod
    def merge_files(merge, paths):
        out = file(merge, 'w')
        for path in paths:
            temp = file(path, 'r')
            for line in temp:
                out.write(line)
        return merge

    @staticmethod
    def merge_fastas(merged_path, fasta_paths):
        for path in fasta_path:
            merge_files(merged_path, fasta_path)

        return merged_path

    @staticmethod
    def assert_file(path, cmd = None):
        if not os.path.exists(path) or os.path.getsize(path) == 0:
            print >> sys.stderr, "###ERROR"
            print >> sys.stderr, "File not created: ", path
            if not cmd == None:
                print >> sys.stderr, "From command: ", cmd
            print >> sys.stderr, "###"
            sys.exit(-1)


class samtools:
    @staticmethod
    def index_fasta(fasta):
        return common.system("samtools faidx {}".format(fasta))

    @staticmethod
    def sam2bam(sam, bam):
        return common.system("samtools view -bS -o {} {}".format(bam, sam))

    @staticmethod
    def sortbam(bam, prefix):
        return common.system("samtools sort {} {}".format(bam, prefix))

    @staticmethod
    def index(bam_path):
        cmd = "samtools index {}".format(bam_path)
        common.system(cmd)
        common.assert_file(bam_path, cmd)
        return bam_path

class picard_tools:
    _PICARD_TOOLS_DIR = "~/local/share/picard_tools/"
    _ADD_OR_REPLACE_READ_GROUPS = os.path.join(_PICARD_TOOLS_DIR, "AddOrReplaceReadGroups.jar")
    _MARK_DUPLICATES            = os.path.join(_PICARD_TOOLS_DIR, "MarkDuplicates.jar") 
    _SORT_SAM                   = os.path.join(_PICARD_TOOLS_DIR, "SortSam.jar") 

    @staticmethod
    def mark_duplicates(bam_path, output_bam_path, output_metrics_path):
        cmd = "java -Xmx1g -jar {} \
              INPUT={} \
              OUTPUT={} \
              METRICS_FILE={}"\
              .format(picard_tools._MARK_DUPLICATES, bam_path, output_bam_path, output_metrics_path)
        common.system(cmd)
        common.assert_file(output_bam_path, cmd)
        common.assert_file(output_metrics_path, cmd)
    
        return output_bam_path
    
    @staticmethod
    def add_or_replace_read_groups(input_path, output_path):
        cmd = "java -jar {} I={} O={}  SORT_ORDER=coordinate LB=FOO PL=ILLUMINA PU=BAR SM=NEE".format(picard_tools._ADD_OR_REPLACE_READ_GROUPS, input_path, output_path)
        common.system(cmd) 
        common.assert_file(output_path, cmd)
        return output_path

    @staticmethod
    def sort_sam(aln_path, output_aln_path, sort_option = "coordinate"):
        cmd = "java -jar {} INPUT={} OUTPUT={} SORT_ORDER={}".format(picard_tools._SORT_SAM, aln_path, output_aln_path, sort_option)
        common.system(cmd)
        common.assert_file(output_aln_path, cmd)
        return output_aln_path


class gatk:
    prefix = "java -Xmx2g -jar {} -T ".format("~/local/share/gatk/GenomeAnalysisTK.jar")

    @staticmethod
    def realigner_target_creator(intervals, fasta, bam):
        cmd = gatk.prefix + "RealignerTargetCreator -R {} -I {} -o {}".format(fasta, bam, intervals)
        common.system(cmd)
        return intervals

    @staticmethod
    def indel_realigner(ref_path, bam_path, realigned_bam_path, intervals_path):
        cmd = gatk.prefix + "IndelRealigner \
              -R {} \
              -I {} \
              -targetIntervals {} \
              -o {}".format(ref_path, bam_path, intervals_path, realigned_bam_path)
        common.system(cmd)
        common.assert_file(realigned_bam_path, cmd)
        return realigned_bam_path

    @staticmethod
    def unified_genotyper(ref_path, bam_path, output_vcf_path):
        cmd = gatk.prefix + "UnifiedGenotyper \
              -R {} \
              -I {} \
              -o {} \
              -l INFO \
              -glm BOTH".format(ref_path, bam_path, output_vcf_path)
        common.system(cmd)
        common.assert_file(output_vcf_path, cmd)
        return output_vcf_path



def do_bwa(ref_paths, read_paths, output_dir, pair_ended):
    output_dir = common.create_path(output_dir)
    done_file = os.path.join(output_dir, "bwa.done")

    if not os.path.exists(done_file):
        #Build new fasta path strings.
        fasta_paths = (os.path.basename(path) for path in ref_paths)
        fasta_paths = (path.replace(".gbk", ".fasta") for path in fasta_paths)
        fasta_paths = [os.path.join(output_dir, path) for path in fasta_paths]
        map(common.gbk2fasta, ref_paths, fasta_paths)

        fasta_path = "merged.fasta" if len(fasta_paths) > 1 else fasta_paths[0]
        common.system("bwa index {}".format(fasta_path))

        #Build new index file paths.
        index_paths = (os.path.basename(path) for path in read_paths)
        index_paths = (path.replace(".fastq", ".sai") for path in index_paths)
        index_paths = [os.path.join(output_dir, path) for path in index_paths]

        for index, read in zip(index_paths, read_paths):
            common.system("bwa aln -f {} {} {}".format(index, fasta_path, read))

        #Create alignment file.
        sam_path = os.path.join(output_dir, "reference.sam")
        if pair_ended == True:
            common.system("bwa sampe -f {} {} {} {}".format(sam_path, fasta_path, " ".join(index_paths), " ".join(read_paths)))
        else:
            common.system("bwa samse -f {} {} {} {}".format(sam_path, fasta_path, " ".join(index_paths), " ".join(read_paths)))


        p.dump((fasta_path, sam_path), open(done_file, 'w'))
    else:
        fasta_path, sam_path = p.load(open(done_file, 'r'))

    print >> sys.stderr, "BWA End."
    print >> sys.stderr, "\tReference file:", fasta_path
    print >> sys.stderr, "\tAlignment file:", sam_path
    return (fasta_path, sam_path)


def do_bowtie(ref_paths, read_paths, output_dir, pair_ended):
    output_dir = common.create_path(output_dir)
    done_file = os.path.join(output_dir, "bowtie.done")

    if not os.path.exists(done_file):
        #Build new fasta path strings.
        fasta_paths = (os.path.basename(path) for path in ref_paths)
        fasta_paths = (path.replace(".gbk", ".fasta") for path in fasta_paths)
        fasta_paths = [os.path.join(output_dir, path) for path in fasta_paths]
        map(common.gbk2fasta, ref_paths, fasta_paths)

        fasta_path = "merged.fasta" if len(fasta_paths) > 1 else fasta_paths[0]

        #Index
        prefix = fasta_path.replace(".fasta", "")
        common.system("bowtie-build {} {}".format(fasta_path, prefix))

        #Create alignment file.
        sam_path = os.path.join(output_dir, "reference.sam")
        fastq_ = ""
        if pair_ended == True:
            """ 
            Requires matched -1 file_1.fastq -2 file_2.fastq parameters.
            Works under the assumption that user passed the read parameters in order.
            """
            fastq_ = " ".join(("-{} {}".format((i % 2) + 1, path) for i, path in enumerate(read_paths)))
        else:
            fastq_ = ",".join(read_paths)

        common.system("bowtie -t --sam {} {} {}".format(prefix, fastq_, sam_path))

        p.dump((fasta_path, sam_path), open(done_file, 'w'))
    else:
        fasta_path, sam_path = p.load(open(done_file, 'r'))

    print >> sys.stderr, "Bowtie End."
    print >> sys.stderr, "\tReference file:", fasta_path
    print >> sys.stderr, "\tAlignment file:", sam_path
    return (fasta_path, sam_path)



def do_gatk(args):
    output_dir = common.create_path(args.output_dir)

    mapper = args.mapper.lower()
    if mapper == "bowtie":
        fasta, sam = do_bowtie(args.ref_paths, args.read_paths, args.output_dir, args.pair_ended)
    elif mapper == "bwa":
        fasta, sam = do_bwa(args.ref_paths, args.read_paths, args.output_dir, args.pair_ended)

    #Add read groups
    RG_sam = os.path.join(output_dir, "RG.sam")
    picard_tools.add_or_replace_read_groups(sam, RG_sam)
    sam = RG_sam

    #Convert to bam
    bam = os.path.join(output_dir, os.path.basename(sam.replace(".sam", ".bam")))
    samtools.sam2bam(sam, bam)

    #MarkDuplicates
    dedup_metrics_path = os.path.join(output_dir, "dedup.metrics")
    dedup_bam = os.path.join(output_dir, "dedup.bam")
    picard_tools.mark_duplicates(bam, dedup_bam, dedup_metrics_path)
    bam = dedup_bam

    #Index
    samtools.index(bam)

    #Realign targets
    intervals = os.path.join(output_dir, "output.intervals")
    gatk.realigner_target_creator(intervals, fasta, bam)

    #Apply realignment
    realigned_bam = os.path.join(output_dir, "realigned.bam")
    gatk.indel_realigner(fasta, bam, realigned_bam, intervals)
    bam = realigned_bam

    #Call variants
    vcf_path = os.path.join(output_dir, "output.vcf")
    gatk.unified_genotyper(fasta, bam, vcf_path)

    print >> sys.stderr, "GATK created : " + vcf_path
    print >> sys.stderr, "Suggested filters for SNPs:"
    print >> sys.stderr, "\tQD < 2.0"
    print >> sys.stderr, "\tMQ < 40.0"
    print >> sys.stderr, "\tFS > 60.0"
    print >> sys.stderr, "\tHaplotypeScore > 13.0"
    print >> sys.stderr, "\tMQRankSum < -12.5"
    print >> sys.stderr, "\tReadPosRankSum < -8.0"
    print >> sys.stderr, ""
    print >> sys.stderr, "Suggested filters for INDELs:"
    print >> sys.stderr, "\tQD < 2.0"
    print >> sys.stderr, "\tReadPosRankSum < -20.0"
    print >> sys.stderr, "\tInbreedingCoeff < -0.8"
    print >> sys.stderr, "\tFS > 200.0"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-r", action = "append", dest = "ref_paths", help = "GBK format", required = True)
    parser.add_argument("-o", dest = "output_dir", default = "gatk")
    parser.add_argument("--pe", action = "store_true", dest = "pair_ended", default = False)
    parser.add_argument("--mapper", dest = "mapper", default = "bwa", choices = ["bwa", "bowtie"])
    parser.add_argument("read_paths", nargs = '+')
    parser.set_defaults(func = do_gatk)

    args =parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
