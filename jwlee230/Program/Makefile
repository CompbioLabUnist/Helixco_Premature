# Variables
DATE := $(shell date "+%y%m%d")
RANDOM := $(shell bash -c 'echo $$RANDOM')
IMAGE_NAME = premature:latest
SHELL := /bin/bash
CPUS = 30
MEMS = 200G
PWD := $(shell pwd)
TOOLS =
DOCKER = docker run --rm $(VOLUME_OPTS) $(RUN_OPTS) $(IMAGE_NAME)

# Options
VOLUME_OPTS = --volume $(abspath Output):/Output --volume $(abspath Data):/Data:ro --volume $(abspath Python):/Python --volume $(abspath Bash):/Bash:ro --tmpfs /tmpfs
RUN_OPTS = --tty --cpus="$(CPUS)" --memory="$(MEMS)"

# General
all:
.PHONY: all

log Output Docker/Tools:
	mkdir $@

# Tools

# Docker
TOOLS += $(wildcard Docker/*)
build.log: Docker/Dockerfile $(TOOLS) | log Output
	rm -fv $@
	docker images | grep $(IMAGE_NAME) && docker rmi $(IMAGE_NAME) || true
	docker build --rm --tag $(IMAGE_NAME) $(<D) | tee $@

build: build.log
.PHONY: build

interactive: build.log
	docker run --rm $(VOLUME_OPTS) $(RUN_OPTS) --interactive $(IMAGE_NAME) /bin/bash || true
.PHONY: interactive

delete: build.log
	docker rmi $(IMAGE_NAME)
	rm -fv build.log
.PHONY: delete

stop:
	docker rm $(CONTAINER_NAME)

# SGE
run: | log Output
	echo -e "#!/bin/bash\nmake -C $(PWD) latest" > tmp.sh
	sbatch --chdir=$(abspath .) --cpus-per-task=$(CPUS) --error="$(abspath log)/%x.e%A" --job-name="Premature_$(DATE)" --mail-type=ALL --mail-user="230@fumire.moe" --mem=$(MEMS) --output="$(abspath log)/%x.o%A" --export=ALL tmp.sh
.PHONY: run

# Actual
latest: step01 step02 step03 step04 step05 step06 step07 step08 step09 step21 step22 step26 step28 step29 step30 step31 step32 step33 step34 step35 step36
.PHONY: latest

# Step 01 (Make manifest files)
Output/Manifest:
	mkdir -p $@

Output/Manifest/first.tsv: Python/step01.py $(wildcard Data/RawData/*.fastq.gz) | Output/Manifest build.log
	$(DOCKER) python3 $(addprefix /,$^) > $@

Output/Manifest/second.tsv: Python/step01.py $(wildcard Data/HN00145920/*.fastq.gz) | Output/Manifest build.log
	$(DOCKER) python3 $(addprefix /,$^) > $@

Output/Manifest/EBI.tsv: Python/step01-2.py $(wildcard Data/EBIData/ena_files/*/*.fastq.gz) | Output/Manifest build.log
	$(DOCKER) python3 $(addprefix /,$^) > $@

Output/Manifest/HMP.tsv: Python/step01-3.py $(wildcard Data/HMPData/*.fastq.gz) | Output/Manifest build.log
	$(DOCKER) python3 $(addprefix /,$^) > $@

Output/Manifest/Stool.tsv: Python/step01-4.py Data/Metadata/readable.xlsx $(wildcard Data/StoolData/*.fastq.gz) | Output/Manifest build.log
	$(DOCKER) python3 $(addprefix /,$(filter %.py,$^) $(filter %.xlsx,$^)) $(foreach f,$(wildcard Data/StoolData/*.fastq.gz),--stool $(addprefix /,$f)) > $@

Output/Manifest/everything.tsv: Python/step01-5.py Data/Metadata/readable.xlsx $(wildcard Data/RawData/*.fastq.gz) $(wildcard Data/HN00145920/*.fastq.gz) $(wildcard Data/HN00148220/*.fastq.gz) $(wildcard Data/StoolData/*.fastq.gz) | Output/Manifest build.log
	$(DOCKER) python3 $(addprefix /,$(filter %.py,$^) $(filter %.xlsx,$^)) $(foreach f,$(wildcard Data/RawData/*.fastq.gz),--first $(addprefix /,$f)) $(foreach f,$(wildcard Data/HN00145920/*.fastq.gz),--second $(addprefix /,$f)) $(foreach f,$(wildcard Data/HN00148220/*.fastq.gz),--third $(addprefix /,$f)) > $@

step01: Output/Manifest/Stool.tsv Output/Manifest/everything.tsv
.PHONY: step01

# Step 02 (De-multiplexing)
Output/Demux:
	mkdir -p $@

Output/Demux/first.qza: Output/Manifest/first.tsv | Output/Demux build.log
	$(DOCKER) qiime tools import --type "SampleData[PairedEndSequencesWithQuality]" --input-format "PairedEndFastqManifestPhred33V2" --input-path $(addprefix /,$<) --output-path $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

Output/Demux/second.qza: Output/Manifest/second.tsv | Output/Demux build.log
	$(DOCKER) qiime tools import --type "SampleData[PairedEndSequencesWithQuality]" --input-format "PairedEndFastqManifestPhred33V2" --input-path $(addprefix /,$<) --output-path $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

Output/Demux/EBI.qza: Output/Manifest/EBI.tsv | Output/Demux build.log
	$(DOCKER) qiime tools import --type "SampleData[SequencesWithQuality]" --input-format "SingleEndFastqManifestPhred33V2" --input-path $(addprefix /,$<) --output-path $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

Output/Demux/HMP.qza: Output/Manifest/HMP.tsv | Output/Demux build.log
	$(DOCKER) qiime tools import --type "SampleData[PairedEndSequencesWithQuality]" --input-format "PairedEndFastqManifestPhred33V2" --input-path $(addprefix /,$<) --output-path $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

Output/Demux/Stool.qza: Output/Manifest/Stool.tsv | Output/Demux build.log
	$(DOCKER) qiime tools import --type "SampleData[PairedEndSequencesWithQuality]" --input-format "PairedEndFastqManifestPhred33V2" --input-path $(addprefix /,$<) --output-path $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

Output/Demux/everything.qza: Output/Manifest/everything.tsv | Output/Demux build.log
	$(DOCKER) qiime tools import --type "SampleData[PairedEndSequencesWithQuality]" --input-format "PairedEndFastqManifestPhred33V2" --input-path $(addprefix /,$<) --output-path $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

Output/Demux/%.qzv: Output/Demux/%.qza | build.log
	$(DOCKER) qiime demux summarize --i-data $(addprefix /,$<) --o-visualization $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

step02: Output/Demux/Stool.qza Output/Demux/Stool.qzv Output/Demux/everything.qza Output/Demux/everything.qzv
.PHONY: step02

# Step 03 (Quality filter)
Output/Filter:
	mkdir -p $@

Output/Filter/%.seq.qza Output/Filter/%.stat.qza: Output/Demux/%.qza | Output/Filter build.log
	$(DOCKER) qiime quality-filter q-score --i-demux $(addprefix /,$<) --o-filtered-sequences $(addprefix /,$(basename $(basename $@))).seq.qza --o-filter-stats $(addprefix /,$(basename $(basename $@))).stat.qza 1> $@.stdout 2> $@.stderr

step03: Output/Filter/Stool.seq.qza Output/Filter/everything.seq.qza
.PHONY: step03

# Step 04 (DADA2 Denoising)
Output/DADA2:
	mkdir -p $@

Output/DADA2/first.table.qza Output/DADA2/first.seq.qza Output/DADA2/first.stat.qza: Output/Demux/first.qza | Output/DADA2 build.log
	@echo "You have to check the length manually!!"
	$(DOCKER) qiime dada2 denoise-paired --i-demultiplexed-seqs $(addprefix /,$<) --p-n-threads $(CPUS) --p-trunc-len-f "300" --p-trunc-len-r "265" --o-table $(addprefix /,$(basename $(basename $@))).table.qza --o-representative-sequences $(addprefix /,$(basename $(basename $@))).seq.qza --o-denoising-stats $(addprefix /,$(basename $(basename $@))).stat.qza 1> $@.stdout 2> $@.stderr

Output/DADA2/second.table.qza Output/DADA2/second.seq.qza Output/DADA2/second.stat.qza: Output/Demux/second.qza | Output/DADA2 build.log
	@echo "You have to check the length manually!!"
	$(DOCKER) qiime dada2 denoise-paired --i-demultiplexed-seqs $(addprefix /,$<) --p-n-threads $(CPUS) --p-trunc-len-f "300" --p-trunc-len-r "222" --o-table $(addprefix /,$(basename $(basename $@))).table.qza --o-representative-sequences $(addprefix /,$(basename $(basename $@))).seq.qza --o-denoising-stats $(addprefix /,$(basename $(basename $@))).stat.qza 1> $@.stdout 2> $@.stderr

Output/DADA2/EBI.table.qza Output/DADA2/EBI.seq.qza Output/DADA2/EBI.stat.qza: Output/Demux/EBI.qza | Output/DADA2 build.log
	@echo "You have to check the length manually!!"
	$(DOCKER) qiime dada2 denoise-single --i-demultiplexed-seqs $(addprefix /,$<) --p-trunc-len "150" --p-n-threads $(CPUS) --o-table $(addprefix /,$(basename $(basename $@))).table.qza --o-representative-sequences $(addprefix /,$(basename $(basename $@))).seq.qza --o-denoising-stats $(addprefix /,$(basename $(basename $@))).stat.qza 1> $@.stdout 2> $@.stderr

Output/DADA2/HMP.table.qza Output/DADA2/HMP.seq.qza Output/DADA2/HMP.stat.qza: Output/Demux/HMP.qza | Output/DADA2 build.log
	@echo "You have to check the length manually!!"
	$(DOCKER) qiime dada2 denoise-paired --i-demultiplexed-seqs $(addprefix /,$<) --p-n-threads $(CPUS) --p-trunc-len-f "278" --p-trunc-len-r "226" --o-table $(addprefix /,$(basename $(basename $@))).table.qza --o-representative-sequences $(addprefix /,$(basename $(basename $@))).seq.qza --o-denoising-stats $(addprefix /,$(basename $(basename $@))).stat.qza 1> $@.stdout 2> $@.stderr

Output/DADA2/Stool.table.qza Output/DADA2/Stool.seq.qza Output/DADA2/Stool.stat.qza: Output/Demux/Stool.qza | Output/DADA2 build.log
	@echo "You have to check the length manually!!"
	$(DOCKER) qiime dada2 denoise-paired --i-demultiplexed-seqs $(addprefix /,$<) --p-n-threads $(CPUS) --p-trunc-len-f "250" --p-trunc-len-r "251" --o-table $(addprefix /,$(basename $(basename $@))).table.qza --o-representative-sequences $(addprefix /,$(basename $(basename $@))).seq.qza --o-denoising-stats $(addprefix /,$(basename $(basename $@))).stat.qza 1> $@.stdout 2> $@.stderr

Output/DADA2/everything.table.qza Output/DADA2/everything.seq.qza Output/DADA2/everything.stat.qza: Output/Demux/everything.qza | Output/DADA2 build.log
	@echo "You have to check the length manually!!"
	$(DOCKER) qiime dada2 denoise-paired --i-demultiplexed-seqs $(addprefix /,$<) --p-n-threads $(CPUS) --p-trunc-len-f "300" --p-trunc-len-r "245" --o-table $(addprefix /,$(basename $(basename $@))).table.qza --o-representative-sequences $(addprefix /,$(basename $(basename $@))).seq.qza --o-denoising-stats $(addprefix /,$(basename $(basename $@))).stat.qza 1> $@.stdout 2> $@.stderr

step04: Output/DADA2/Stool.seq.qza Output/DADA2/everything.seq.qza
.PHONY: step04

# Step 05 (Deblur Denoising)
Output/Deblur:
	mkdir -p $@

Output/Deblur/first.table.qza Output/Deblur/first.seq.qza Output/Deblur/first.stat.qza: Output/Filter/first.seq.qza | Output/Deblur build.log
	@echo "You have to check the length manually!!"
	$(DOCKER) qiime deblur denoise-16S --i-demultiplexed-seqs $(addprefix /,$<) --p-trim-length "265" --p-jobs-to-start $(CPUS) --o-table $(addprefix /,$(basename $(basename $@))).table.qza --o-representative-sequences $(addprefix /,$(basename $(basename $@))).seq.qza --o-stats $(addprefix /,$(basename $(basename $@))).stat.qza 1> $@.stdout 2> $@.stderr

Output/Deblur/second.table.qza Output/Deblur/second.seq.qza Output/Deblur/second.stat.qza: Output/Filter/second.seq.qza | Output/Deblur build.log
	@echo "You have to check the length manually!!"
	$(DOCKER) qiime deblur denoise-16S --i-demultiplexed-seqs $(addprefix /,$<) --p-trim-length "222" --p-jobs-to-start $(CPUS) --o-table $(addprefix /,$(basename $(basename $@))).table.qza --o-representative-sequences $(addprefix /,$(basename $(basename $@))).seq.qza --o-stats $(addprefix /,$(basename $(basename $@))).stat.qza 1> $@.stdout 2> $@.stderr

Output/Deblur/EBI.table.qza Output/Deblur/EBI.seq.qza Output/Deblur/EBI.stat.qza: Output/Filter/EBI.seq.qza | Output/Deblur build.log
	@echo "You have to check the length manually!!"
	$(DOCKER) qiime deblur denoise-16S --i-demultiplexed-seqs $(addprefix /,$<) --p-trim-length "150" --p-jobs-to-start $(CPUS) --o-table $(addprefix /,$(basename $(basename $@))).table.qza --o-representative-sequences $(addprefix /,$(basename $(basename $@))).seq.qza --o-stats $(addprefix /,$(basename $(basename $@))).stat.qza 1> $@.stdout 2> $@.stderr

Output/Deblur/HMP.table.qza Output/Deblur/HMP.seq.qza Output/Deblur/HMP.stat.qza: Output/Filter/HMP.seq.qza | Output/Deblur build.log
	@echo "You have to check the length manually!!"
	$(DOCKER) qiime deblur denoise-16S --i-demultiplexed-seqs $(addprefix /,$<) --p-trim-length "226" --p-jobs-to-start $(CPUS) --o-table $(addprefix /,$(basename $(basename $@))).table.qza --o-representative-sequences $(addprefix /,$(basename $(basename $@))).seq.qza --o-stats $(addprefix /,$(basename $(basename $@))).stat.qza 1> $@.stdout 2> $@.stderr

Output/Deblur/Stool.table.qza Output/Deblur/Stool.seq.qza Output/Deblur/Stool.stat.qza: Output/Filter/Stool.seq.qza | Output/Deblur build.log
	@echo "You have to check the length manually!!"
	$(DOCKER) qiime deblur denoise-16S --i-demultiplexed-seqs $(addprefix /,$<) --p-trim-length "250" --p-jobs-to-start $(CPUS) --o-table $(addprefix /,$(basename $(basename $@))).table.qza --o-representative-sequences $(addprefix /,$(basename $(basename $@))).seq.qza --o-stats $(addprefix /,$(basename $(basename $@))).stat.qza 1> $@.stdout 2> $@.stderr

Output/Deblur/everything.table.qza Output/Deblur/everything.seq.qza Output/Deblur/everything.stat.qza: Output/Filter/everything.seq.qza | Output/Deblur build.log
	@echo "You have to check the length manually!!"
	$(DOCKER) qiime deblur denoise-16S --i-demultiplexed-seqs $(addprefix /,$<) --p-trim-length "245" --p-jobs-to-start $(CPUS) --o-table $(addprefix /,$(basename $(basename $@))).table.qza --o-representative-sequences $(addprefix /,$(basename $(basename $@))).seq.qza --o-stats $(addprefix /,$(basename $(basename $@))).stat.qza 1> $@.stdout 2> $@.stderr

step05: Output/Deblur/Stool.seq.qza Output/Deblur/everything.seq.qza
.PHONY: step05

# Step 06 (Taxonomy information)
Output/Taxonomy:
	mkdir -p $@

Output/Taxonomy/silva.qza: | Output/Taxonomy build.log
	wget "https://data.qiime2.org/2021.4/common/silva-138-99-nb-classifier.qza" -O $@
	touch $@

Output/Taxonomy/gg.qza: | Output/Taxonomy build.log
	wget "https://data.qiime2.org/2021.4/common/gg-13-8-99-nb-classifier.qza" -O $@
	touch $@

Output/Taxonomy/homd.fasta: | Output/Taxonomy
	wget "http://www.homd.org/ftp/16S_rRNA_refseq/HOMD_16S_rRNA_RefSeq/V15.22/HOMD_16S_rRNA_RefSeq_V15.22.fasta" -O $@
	touch $@

Output/Taxonomy/homd.txt: | Output/Taxonomy
	wget "http://www.homd.org/ftp/16S_rRNA_refseq/HOMD_16S_rRNA_RefSeq/V15.22/HOMD_16S_rRNA_RefSeq_V15.22.qiime.taxonomy" -O $@
	touch $@

Output/Taxonomy/homd_fasta.qza: Output/Taxonomy/homd.fasta | build.log
	$(DOCKER) qiime tools import --type "FeatureData[Sequence]" --input-path $(addprefix /,$<) --output-path $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

Output/Taxonomy/homd_taxo.qza: Output/Taxonomy/homd.txt | build.log
	$(DOCKER) qiime tools import --type "FeatureData[Taxonomy]" --input-format "HeaderlessTSVTaxonomyFormat" --input-path $(addprefix /,$<) --output-path $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

Output/Taxonomy/homd_seqs.qza: Output/Taxonomy/homd_fasta.qza | build.log
	$(DOCKER) qiime feature-classifier extract-reads --i-sequences $(addprefix /,$<) --p-f-primer "GTGCCAGCMGCCGCGGTAA" --p-r-primer "GGACTACHVGGGTWTCTAAT" --o-reads $(addprefix /,$@) --p-n-jobs $(CPUS) 1> $@.stdout 2> $@.stderr

Output/Taxonomy/homd.qza: Output/Taxonomy/homd_seqs.qza Output/Taxonomy/homd_taxo.qza | build.log
	$(DOCKER) qiime feature-classifier fit-classifier-naive-bayes --i-reference-reads $(addprefix /,$(word 1,$^)) --i-reference-taxonomy $(addprefix /,$(word 2,$^)) --o-classifier $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

Output/Taxonomy/%.DADA2.silva.qza: Output/Taxonomy/silva.qza Output/DADA2/%.seq.qza | build.log
	$(DOCKER) qiime feature-classifier classify-sklearn --i-classifier $(addprefix /,$(word 1,$^)) --i-reads $(addprefix /,$(word 2,$^)) --p-n-jobs $(CPUS) --o-classification $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

Output/Taxonomy/%.DADA2.gg.qza: Output/Taxonomy/gg.qza Output/DADA2/%.seq.qza | build.log
	$(DOCKER) qiime feature-classifier classify-sklearn --i-classifier $(addprefix /,$(word 1,$^)) --i-reads $(addprefix /,$(word 2,$^)) --p-n-jobs $(CPUS) --o-classification $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

Output/Taxonomy/%.DADA2.homd.qza: Output/Taxonomy/homd.qza Output/DADA2/%.seq.qza | build.log
	$(DOCKER) qiime feature-classifier classify-sklearn --i-classifier $(addprefix /,$(word 1,$^)) --i-reads $(addprefix /,$(word 2,$^)) --p-n-jobs $(CPUS) --o-classification $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

Output/Taxonomy/%.Deblur.silva.qza: Output/Taxonomy/silva.qza Output/Deblur/%.seq.qza | build.log
	$(DOCKER) qiime feature-classifier classify-sklearn --i-classifier $(addprefix /,$(word 1,$^)) --i-reads $(addprefix /,$(word 2,$^)) --p-n-jobs $(CPUS) --o-classification $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

Output/Taxonomy/%.Deblur.gg.qza: Output/Taxonomy/gg.qza Output/Deblur/%.seq.qza | build.log
	$(DOCKER) qiime feature-classifier classify-sklearn --i-classifier $(addprefix /,$(word 1,$^)) --i-reads $(addprefix /,$(word 2,$^)) --p-n-jobs $(CPUS) --o-classification $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

Output/Taxonomy/%.Deblur.homd.qza: Output/Taxonomy/homd.qza Output/Deblur/%.seq.qza | build.log
	$(DOCKER) qiime feature-classifier classify-sklearn --i-classifier $(addprefix /,$(word 1,$^)) --i-reads $(addprefix /,$(word 2,$^)) --p-n-jobs $(CPUS) --o-classification $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

step06-reference: Output/Taxonomy/silva.qza Output/Taxonomy/gg.qza Output/Taxonomy/homd.fasta Output/Taxonomy/homd.txt Output/Taxonomy/homd_fasta.qza Output/Taxonomy/homd_taxo.qza Output/Taxonomy/homd_seqs.qza Output/Taxonomy/homd.qza
step06-mapping: Output/Taxonomy/everything.DADA2.gg.qza Output/Taxonomy/everything.DADA2.silva.qza Output/Taxonomy/everything.DADA2.homd.qza Output/Taxonomy/everything.Deblur.gg.qza Output/Taxonomy/everything.Deblur.silva.qza Output/Taxonomy/everything.Deblur.homd.qza Output/Taxonomy/Stool.DADA2.gg.qza Output/Taxonomy/Stool.DADA2.silva.qza Output/Taxonomy/Stool.DADA2.homd.qza Output/Taxonomy/Stool.Deblur.gg.qza Output/Taxonomy/Stool.Deblur.silva.qza Output/Taxonomy/Stool.Deblur.homd.qza
step06: step06-reference step06-mapping
.PHONY: step06 step06-reference step06-mapping

# Step 07 (Build Metadata)
Output/Metadata:
	mkdir -p $@

Output/Metadata/everything.tsv: Python/step07-4.py Output/Manifest/everything.tsv Data/Metadata/metadata.xlsx | Output/Metadata build.log
	$(DOCKER) python3 $(addprefix /,$^) > $@

Output/Metadata/Stool.tsv: Python/step07-4.py Output/Manifest/Stool.tsv Data/Metadata/metadata.xlsx | Output/Metadata build.log
	$(DOCKER) python3 $(addprefix /,$^) > $@

step07: Output/Metadata/everything.tsv Output/Metadata/Stool.tsv
.PHONY: step07

# Step 08 (Run ANCOM)
Output/ANCOM:
	mkdir -p $@

Output/ANCOM/%.DADA2.gg.qza: Output/DADA2/%.table.qza Output/Taxonomy/%.DADA2.gg.qza | Output/ANCOM build.log
	$(DOCKER) qiime taxa collapse --i-table $(addprefix /,$(word 1,$^)) --i-taxonomy $(addprefix /,$(word 2,$^)) --p-level 7 --o-collapsed-table $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

Output/ANCOM/%.DADA2.silva.qza: Output/DADA2/%.table.qza Output/Taxonomy/%.DADA2.silva.qza | Output/ANCOM build.log
	$(DOCKER) qiime taxa collapse --i-table $(addprefix /,$(word 1,$^)) --i-taxonomy $(addprefix /,$(word 2,$^)) --p-level 7 --o-collapsed-table $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

Output/ANCOM/%.Deblur.gg.qza: Output/Deblur/%.table.qza Output/Taxonomy/%.Deblur.gg.qza | Output/ANCOM build.log
	$(DOCKER) qiime taxa collapse --i-table $(addprefix /,$(word 1,$^)) --i-taxonomy $(addprefix /,$(word 2,$^)) --p-level 7 --o-collapsed-table $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

Output/ANCOM/%.Deblur.silva.qza: Output/Deblur/%.table.qza Output/Taxonomy/%.Deblur.silva.qza | Output/ANCOM build.log
	$(DOCKER) qiime taxa collapse --i-table $(addprefix /,$(word 1,$^)) --i-taxonomy $(addprefix /,$(word 2,$^)) --p-level 7 --o-collapsed-table $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

Output/ANCOM/everything.%.Mouth.qza: Output/ANCOM/everything.%.qza Output/Metadata/everything.tsv | build.log
	$(DOCKER) qiime feature-table filter-samples --i-table $(addprefix /,$(word 1,$^)) --m-metadata-file $(addprefix /,$(word 2,$^)) --p-where "[Site]='Mouth'" --o-filtered-table $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

Output/ANCOM/everything.%.Cervix.qza: Output/ANCOM/everything.%.qza Output/Metadata/everything.tsv | build.log
	$(DOCKER) qiime feature-table filter-samples --i-table $(addprefix /,$(word 1,$^)) --m-metadata-file $(addprefix /,$(word 2,$^)) --p-where "[Site]='Cervix'" --o-filtered-table $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

Output/ANCOM/everything.%.Vagina.qza: Output/ANCOM/everything.%.qza Output/Metadata/everything.tsv | build.log
	$(DOCKER) qiime feature-table filter-samples --i-table $(addprefix /,$(word 1,$^)) --m-metadata-file $(addprefix /,$(word 2,$^)) --p-where "[Site]='Vagina'" --o-filtered-table $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

Output/ANCOM/everything.%.Neonate.qza: Output/ANCOM/everything.%.qza Output/Metadata/everything.tsv | build.log
	$(DOCKER) qiime feature-table filter-samples --i-table $(addprefix /,$(word 1,$^)) --m-metadata-file $(addprefix /,$(word 2,$^)) --p-where "[Site] IN ('Neonate-1day', 'Neonate-3day', 'Neonate-5day')" --o-filtered-table $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

Output/ANCOM/everything.%.Neonate-1day.qza: Output/ANCOM/everything.%.qza Output/Metadata/everything.tsv | build.log
	$(DOCKER) qiime feature-table filter-samples --i-table $(addprefix /,$(word 1,$^)) --m-metadata-file $(addprefix /,$(word 2,$^)) --p-where "[Site]=='Neonate-1day'" --o-filtered-table $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

Output/ANCOM/everything.%.Neonate-3day.qza: Output/ANCOM/everything.%.qza Output/Metadata/everything.tsv | build.log
	$(DOCKER) qiime feature-table filter-samples --i-table $(addprefix /,$(word 1,$^)) --m-metadata-file $(addprefix /,$(word 2,$^)) --p-where "[Site]=='Neonate-3day'" --o-filtered-table $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

Output/ANCOM/everything.%.Neonate-5day.qza: Output/ANCOM/everything.%.qza Output/Metadata/everything.tsv | build.log
	$(DOCKER) qiime feature-table filter-samples --i-table $(addprefix /,$(word 1,$^)) --m-metadata-file $(addprefix /,$(word 2,$^)) --p-where "[Site]=='Neonate-5day'" --o-filtered-table $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

Output/ANCOM/%.pseudocount.qza: Output/ANCOM/%.qza | build.log
	$(DOCKER) qiime composition add-pseudocount --i-table $(addprefix /,$(word 1,$^)) --o-composition-table $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

Output/ANCOM/everything.%.PTB.qzv: Output/ANCOM/everything.%.qza Output/Metadata/everything.tsv | build.log
	$(DOCKER) qiime composition ancom --i-table $(addprefix /,$(word 1,$^)) --m-metadata-file $(addprefix /,$(word 2,$^)) --m-metadata-column "Premature" --o-visualization $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

Output/ANCOM/everything.%.Detail_PTB.qzv: Output/ANCOM/everything.%.qza Output/Metadata/everything.tsv | build.log
	$(DOCKER) qiime composition ancom --i-table $(addprefix /,$(word 1,$^)) --m-metadata-file $(addprefix /,$(word 2,$^)) --m-metadata-column "Detail Premature" --o-visualization $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

Output/ANCOM/%/ancom.tsv Output/ANCOM/%/data.tsv: Output/ANCOM/%.qzv | build.log
	$(DOCKER) rm -rf $(addprefix /,$(@D))
	$(DOCKER) qiime tools export --input-path $(addprefix /,$(word 1,$^)) --output-path $(addprefix /,$(@D)) 1> $(@D).stdout 2> $(@D).stderr

Output/ANCOM/%.pdf: Python/step08.py Output/ANCOM/%/ancom.tsv Output/ANCOM/%/data.tsv | build.log
	$(DOCKER) python3 $(addprefix /,$^ $@) 1> $@.stdout 2> $@.stderr

step08-DADA2-gg: Output/ANCOM/everything.DADA2.gg.Neonate.pseudocount.Detail_PTB.pdf Output/ANCOM/everything.DADA2.gg.Neonate-1day.pseudocount.Detail_PTB.pdf Output/ANCOM/everything.DADA2.gg.Neonate-3day.pseudocount.Detail_PTB.pdf Output/ANCOM/everything.DADA2.gg.Neonate-5day.pseudocount.Detail_PTB.pdf Output/ANCOM/everything.DADA2.gg.Neonate.pseudocount.Detail_PTB/ancom.tsv Output/ANCOM/everything.DADA2.gg.Neonate-3day.pseudocount.Detail_PTB.qzv Output/ANCOM/everything.DADA2.gg.Neonate-5day.pseudocount.qza Output/ANCOM/everything.DADA2.gg.Neonate.pseudocount.qza Output/ANCOM/everything.DADA2.gg.Neonate-5day.pseudocount.Detail_PTB.qzv Output/ANCOM/everything.DADA2.gg.Neonate.qza Output/ANCOM/everything.DADA2.gg.Neonate-3day.qza Output/ANCOM/everything.DADA2.gg.Neonate-1day.pseudocount.Detail_PTB.qzv Output/ANCOM/everything.DADA2.gg.Neonate-1day.pseudocount.Detail_PTB/ancom.tsv Output/ANCOM/everything.DADA2.gg.Neonate-3day.pseudocount.Detail_PTB/ancom.tsv Output/ANCOM/everything.DADA2.gg.Neonate-5day.qza Output/ANCOM/everything.DADA2.gg.Neonate.pseudocount.Detail_PTB.qzv Output/ANCOM/everything.DADA2.gg.Neonate-3day.pseudocount.qza Output/ANCOM/everything.DADA2.gg.Neonate-1day.pseudocount.qza Output/ANCOM/everything.DADA2.gg.Neonate-5day.pseudocount.Detail_PTB/ancom.tsv Output/ANCOM/everything.DADA2.gg.Neonate-1day.qza Output/ANCOM/everything.DADA2.gg.Neonate.pseudocount.Detail_PTB/data.tsv Output/ANCOM/everything.DADA2.gg.Neonate-1day.pseudocount.Detail_PTB/data.tsv Output/ANCOM/everything.DADA2.gg.Neonate-3day.pseudocount.Detail_PTB/data.tsv Output/ANCOM/everything.DADA2.gg.Neonate-5day.pseudocount.Detail_PTB/data.tsv
step08-DADA2-silva:
step08-Deblur-gg:
step08-Deblur-silva:
step08: step08-DADA2-gg step08-DADA2-silva step08-Deblur-gg step08-Deblur-silva
.PHONY: step08 step08-DADA2-gg step08-DADA2-silva step08-Deblur-gg step08-Deblur-silva

# Step 09 (Get TSV)
Output/TSV:
	mkdir -p $@

Output/TSV/Stool.DADA2.%.tsv: Bash/export_TSV.bash Output/DADA2/Stool.table.qza Output/Taxonomy/Stool.DADA2.%.qza | Output/TSV build.log
	$(DOCKER) /bin/bash $(addprefix /,$^ $@) 1> $@.stdout 2> $@.stderr

Output/TSV/Stool.Deblur.%.tsv: Bash/export_TSV.bash Output/Deblur/Stool.table.qza Output/Taxonomy/Stool.Deblur.%.qza | Output/TSV build.log
	$(DOCKER) /bin/bash $(addprefix /,$^ $@) 1> $@.stdout 2> $@.stderr

Output/TSV/everything.DADA2.%.tsv: Bash/export_TSV.bash Output/DADA2/everything.table.qza Output/Taxonomy/everything.DADA2.%.qza | Output/TSV build.log
	$(DOCKER) /bin/bash $(addprefix /,$^ $@) 1> $@.stdout 2> $@.stderr

Output/TSV/everything.Deblur.%.tsv: Bash/export_TSV.bash Output/Deblur/everything.table.qza Output/Taxonomy/everything.Deblur.%.qza | Output/TSV build.log
	$(DOCKER) /bin/bash $(addprefix /,$^ $@) 1> $@.stdout 2> $@.stderr

Output/TSV/Macrogen.tsv: Python/step09.py Data/Macrogen/First/otu_table_shared_with_tax_assignment.xlsx Data/Macrogen/Second+Third/Mouth/otu_table_shared_with_tax_assignment.xlsx Data/Macrogen/Second+Third/Vagina/otu_table_shared_with_tax_assignment.xlsx Data/Macrogen/Second+Third/Neonate-1day/otu_table_shared_with_tax_assignment.xlsx Data/Macrogen/Second+Third/Neonate-3day/otu_table_shared_with_tax_assignment.xlsx Data/Macrogen/Second+Third/Neonate-5day/otu_table_shared_with_tax_assignment.xlsx Data/Metadata/readable.xlsx | Output/TSV build.log
	$(DOCKER) python3 $(addprefix /,$^ $@) 1> $@.stdout 2> $@.stderr

step09: Output/TSV/everything.DADA2.gg.tsv Output/TSV/everything.DADA2.silva.tsv Output/TSV/everything.DADA2.homd.tsv Output/TSV/everything.Deblur.gg.tsv Output/TSV/everything.Deblur.silva.tsv Output/TSV/everything.Deblur.homd.tsv Output/TSV/Stool.DADA2.gg.tsv Output/TSV/Stool.DADA2.silva.tsv Output/TSV/Stool.DADA2.homd.tsv Output/TSV/Stool.Deblur.gg.tsv Output/TSV/Stool.Deblur.silva.tsv Output/TSV/Stool.Deblur.homd.tsv Output/TSV/Macrogen.tsv
.PHONY: step09

# Step 21 (Read & clearify first TSV into pandas)
Output/TSV/%.tar.gz: Python/step21.py Output/TSV/%.tsv | build.log
	$(DOCKER) python3 $(addprefix /,$^ $@) 1> $@.stdout 2> $@.stderr

step21: Output/TSV/everything.DADA2.gg.tar.gz Output/TSV/everything.DADA2.silva.tar.gz Output/TSV/everything.DADA2.homd.tar.gz Output/TSV/everything.Deblur.gg.tar.gz Output/TSV/everything.Deblur.silva.tar.gz Output/TSV/everything.Deblur.homd.tar.gz Output/TSV/Stool.DADA2.gg.tar.gz Output/TSV/Stool.DADA2.silva.tar.gz Output/TSV/Stool.DADA2.homd.tar.gz Output/TSV/Stool.Deblur.gg.tar.gz Output/TSV/Stool.Deblur.silva.tar.gz Output/TSV/Stool.Deblur.homd.tar.gz Output/TSV/Macrogen.tar.gz
.PHONY: step21

# Step 22 (Making & Draw t-SNE)
Output/Step22:
	mkdir -p $@

Output/Step22/everything.%.tar.gz: Python/step11.py Output/TSV/everything.%.tsv Output/Metadata/everything.tsv | Output/Step22 build.log
	$(DOCKER) python3 $(addprefix /,$^ $@) --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/Step22/Stool.%.tar.gz: Python/step11.py Output/TSV/Stool.%.tsv Output/Metadata/Stool.tsv | Output/Step22 build.log
	$(DOCKER) python3 $(addprefix /,$^ $@) --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/Step22/Macrogen.tar.gz: Python/step11-2.py Output/TSV/Macrogen.tsv Output/Metadata/everything.tsv | Output/Step22 build.log
	$(DOCKER) python3 $(addprefix /,$^ $@) --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/Step22/%.tar: Python/step22-3.py Output/Step22/%.tar.gz | Output/Step22 build.log
	$(DOCKER) python3 $(addprefix /,$^ $@) 1> $@.stdout 2> $@.stderr

step22: Output/Step22/everything.DADA2.gg.tar.gz Output/Step22/everything.DADA2.gg.tar Output/Step22/everything.DADA2.silva.tar.gz Output/Step22/everything.DADA2.silva.tar Output/Step22/everything.DADA2.homd.tar.gz Output/Step22/everything.DADA2.homd.tar Output/Step22/everything.Deblur.gg.tar.gz Output/Step22/everything.Deblur.gg.tar Output/Step22/everything.Deblur.silva.tar.gz Output/Step22/everything.Deblur.silva.tar Output/Step22/everything.Deblur.homd.tar.gz Output/Step22/everything.Deblur.homd.tar Output/Step22/Stool.DADA2.gg.tar.gz Output/Step22/Stool.DADA2.gg.tar Output/Step22/Stool.DADA2.silva.tar.gz Output/Step22/Stool.DADA2.silva.tar Output/Step22/Stool.Deblur.gg.tar.gz Output/Step22/Stool.Deblur.gg.tar Output/Step22/Stool.Deblur.silva.tar.gz Output/Step22/Stool.Deblur.gg.tar Output/Step22/Macrogen.tar.gz Output/Step22/Macrogen.tar
.PHONY: step22

# Step 26 (RandomForest Classifier)
Output/Step26:
	mkdir -p $@

Output/Step26/RF.%.tar: Python/step26.py Output/TSV/everything.%.tar.gz Output/Metadata/everything.tsv | Output/Step26 build.log
	$(DOCKER) python3 $(addprefix /,$^ $@) --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

step26: Output/Step26/RF.DADA2.gg.tar Output/Step26/RF.DADA2.silva.tar Output/Step26/RF.DADA2.homd.tar Output/Step26/RF.Deblur.gg.tar Output/Step26/RF.Deblur.silva.tar Output/Step26/RF.Deblur.homd.tar
.PHONY: step26

# Step 27 (DecisionTree Classifier)
Output/Step27:
	mkdir -p $@

Output/Step27/DT.%.tar: Python/step27.py Output/TSV/helixco.%.tar.gz Output/TSV/EBI.%.tar.gz Output/TSV/HMP.%.tar.gz Output/Metadata/first.tsv | Output/Step27 build.log
	$(DOCKER) python3 $(addprefix /,$^ $@) --cpu $(CPUS) 1> $@.stdout 2> $@.stderr

step27: Output/Step27/DT.DADA2.gg.tar Output/Step27/DT.DADA2.silva.tar Output/Step27/DT.Deblur.gg.tar Output/Step27/DT.Deblur.silva.tar
.PHONY: step27

# Step 28 (Phylogenetic diversity analysis)
Output/Step28:
	mkdir -p $@

Output/Step28/%.DADA2.aligned-rep-seqs.qza Output/Step28/%.DADA2.masked-aligned-rep-seqs.qza Output/Step28/%.DADA2.unrooted-tree.qza Output/Step28/%.DADA2.rooted-tree.qza: Output/DADA2/%.seq.qza | Output/Step28 build.log
	$(DOCKER) qiime phylogeny align-to-tree-mafft-fasttree --i-sequences $(addprefix /,$<) --o-alignment $(addprefix /,$(basename $(basename $@))).aligned-rep-seqs.qza --o-masked-alignment $(addprefix /,$(basename $(basename $@))).masked-aligned-rep-seqs.qza --o-tree $(addprefix /,$(basename $(basename $@))).unrooted-tree.qza --o-rooted-tree $(addprefix /,$(basename $(basename $@))).rooted-tree.qza --p-n-threads $(CPUS) 1> $@.stdout 2> $@.stderr

Output/Step28/%.Deblur.aligned-rep-seqs.qza Output/Step28/%.Deblur.masked-aligned-rep-seqs.qza Output/Step28/%.Deblur.unrooted-tree.qza Output/Step28/%.Deblur.rooted-tree.qza: Output/Deblur/%.seq.qza | Output/Step28 build.log
	$(DOCKER) qiime phylogeny align-to-tree-mafft-fasttree --i-sequences $(addprefix /,$<) --o-alignment $(addprefix /,$(basename $(basename $@))).aligned-rep-seqs.qza --o-masked-alignment $(addprefix /,$(basename $(basename $@))).masked-aligned-rep-seqs.qza --o-tree $(addprefix /,$(basename $(basename $@))).unrooted-tree.qza --o-rooted-tree $(addprefix /,$(basename $(basename $@))).rooted-tree.qza --p-n-threads $(CPUS) 1> $@.stdout 2> $@.stderr

step28: Output/Step28/everything.DADA2.rooted-tree.qza
.PHONY: step28

# Step 29 (Core phylogenetic metrics)
Output/Step29 Output/Step29/DADA2 Output/Step29/Deblur:
	mkdir -p $@

Output/Step29/%.DADA2.table.qzv: Output/DADA2/%.table.qza Output/Metadata/%.tsv | Output/Step29 build.log
	$(DOCKER) qiime feature-table summarize --i-table $(addprefix /,$(word 1,$^)) --m-sample-metadata-file $(addprefix /,$(word 2,$^)) --o-visualization $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

Output/Step29/%.DADA2.seq.qzv: Output/DADA2/%.seq.qza | Output/Step29 build.log
	$(DOCKER) qiime feature-table tabulate-seqs --i-data $(addprefix /,$<) --o-visualization $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

Output/Step29/%.Deblur.table.qzv: Output/Deblur/%.table.qza Output/Metadata/%.tsv | Output/Step29 build.log
	$(DOCKER) qiime feature-table summarize --i-table $(addprefix /,$(word 1,$^)) --m-sample-metadata-file $(addprefix /,$(word 2,$^)) --o-visualization $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

Output/Step29/%.Deblur.seq.qzv: Output/Deblur/%.seq.qza | Output/Step29 build.log
	$(DOCKER) qiime feature-table tabulate-seqs --i-data $(addprefix /,$<) --o-visualization $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

Output/Step29/DADA2/everything.faith_pd_vector.qza Output/Step29/DADA2/everything.unweighted_unifrac_distance_matrix.qza Output/Step29/DADA2/everything.bray_curtis_pcoa_results.qza Output/Step29/DADA2/everything.shannon_vector.qza Output/Step29/DADA2/everything.rarefied_table.qza Output/Step29/DADA2/everything.weighted_unifrac_distance_matrix.qza Output/Step29/DADA2/everything.jaccard_pcoa_results.qza Output/Step29/DADA2/everything.weighted_unifrac_pcoa_results.qza Output/Step29/DADA2/everything.observed_features_vector.qza Output/Step29/DADA2/everything.jaccard_distance_matrix.qza Output/Step29/DADA2/everything.evenness_vector.qza Output/Step29/DADA2/everything.bray_curtis_distance_matrix.qza Output/Step29/DADA2/everything.unweighted_unifrac_pcoa_results.qza Output/Step29/DADA2/everything.unweighted_unifrac_emperor.qzv Output/Step29/DADA2/everything.jaccard_emperor.qzv Output/Step29/DADA2/everything.bray_curtis_emperor.qzv Output/Step29/DADA2/everything.weighted_unifrac_emperor.qzv: Output/Step28/everything.DADA2.rooted-tree.qza Output/DADA2/everything.table.qza Output/Metadata/everything.tsv | Output/Step29/DADA2 build.log
	@echo "You have to check the sequencing depth manually!!"
	$(DOCKER) qiime diversity core-metrics-phylogenetic --i-phylogeny $(addprefix /,$(word 1,$^)) --i-table $(addprefix /,$(word 2,$^)) --m-metadata-file $(addprefix /,$(word 3,$^)) --p-sampling-depth "1046" --o-rarefied-table $(addprefix /,$(basename $(basename $@))).rarefied_table.qza --o-faith-pd-vector $(addprefix /,$(basename $(basename $@))).faith_pd_vector.qza --o-observed-features-vector $(addprefix /,$(basename $(basename $@))).observed_features_vector.qza --o-shannon-vector $(addprefix /,$(basename $(basename $@))).shannon_vector.qza --o-evenness-vector $(addprefix /,$(basename $(basename $@))).evenness_vector.qza --o-unweighted-unifrac-distance-matrix $(addprefix /,$(basename $(basename $@))).unweighted_unifrac_distance_matrix.qza --o-weighted-unifrac-distance-matrix $(addprefix /,$(basename $(basename $@))).weighted_unifrac_distance_matrix.qza --o-jaccard-distance-matrix $(addprefix /,$(basename $(basename $@))).jaccard_distance_matrix.qza --o-bray-curtis-distance-matrix $(addprefix /,$(basename $(basename $@))).bray_curtis_distance_matrix.qza --o-unweighted-unifrac-pcoa-results $(addprefix /,$(basename $(basename $@))).unweighted_unifrac_pcoa_results.qza --o-weighted-unifrac-pcoa-results $(addprefix /,$(basename $(basename $@))).weighted_unifrac_pcoa_results.qza --o-unweighted-unifrac-emperor $(addprefix /,$(basename $(basename $@))).unweighted_unifrac_emperor.qzv --o-weighted-unifrac-emperor $(addprefix /,$(basename $(basename $@))).weighted_unifrac_emperor.qzv --o-jaccard-emperor $(addprefix /,$(basename $(basename $@))).jaccard_emperor.qzv --o-bray-curtis-emperor $(addprefix /,$(basename $(basename $@))).bray_curtis_emperor.qzv --o-jaccard-pcoa-results $(addprefix /,$(basename $(basename $@))).jaccard_pcoa_results.qza --o-bray-curtis-pcoa-results $(addprefix /,$(basename $(basename $@))).bray_curtis_pcoa_results.qza --p-n-jobs-or-threads $(CPUS) 1> $@.stdout 2> $@.stderr

Output/Step29/Deblur/everything.faith_pd_vector.qza Output/Step29/Deblur/everything.unweighted_unifrac_distance_matrix.qza Output/Step29/Deblur/everything.bray_curtis_pcoa_results.qza Output/Step29/Deblur/everything.shannon_vector.qza Output/Step29/Deblur/everything.rarefied_table.qza Output/Step29/Deblur/everything.weighted_unifrac_distance_matrix.qza Output/Step29/Deblur/everything.jaccard_pcoa_results.qza Output/Step29/Deblur/everything.weighted_unifrac_pcoa_results.qza Output/Step29/Deblur/everything.observed_features_vector.qza Output/Step29/Deblur/everything.jaccard_distance_matrix.qza Output/Step29/Deblur/everything.evenness_vector.qza Output/Step29/Deblur/everything.bray_curtis_distance_matrix.qza Output/Step29/Deblur/everything.unweighted_unifrac_pcoa_results.qza Output/Step29/Deblur/everything.unweighted_unifrac_emperor.qzv Output/Step29/Deblur/everything.jaccard_emperor.qzv Output/Step29/Deblur/everything.bray_curtis_emperor.qzv Output/Step29/Deblur/everything.weighted_unifrac_emperor.qzv: Output/Step28/everything.Deblur.rooted-tree.qza Output/Deblur/everything.table.qza Output/Metadata/everything.tsv | Output/Step29/Deblur build.log
	@echo "You have to check the sequencing depth manually!!"
	$(DOCKER) qiime diversity core-metrics-phylogenetic --i-phylogeny $(addprefix /,$(word 1,$^)) --i-table $(addprefix /,$(word 2,$^)) --m-metadata-file $(addprefix /,$(word 3,$^)) --p-sampling-depth "4864" --o-rarefied-table $(addprefix /,$(basename $(basename $@))).rarefied_table.qza --o-faith-pd-vector $(addprefix /,$(basename $(basename $@))).faith_pd_vector.qza --o-observed-features-vector $(addprefix /,$(basename $(basename $@))).observed_features_vector.qza --o-shannon-vector $(addprefix /,$(basename $(basename $@))).shannon_vector.qza --o-evenness-vector $(addprefix /,$(basename $(basename $@))).evenness_vector.qza --o-unweighted-unifrac-distance-matrix $(addprefix /,$(basename $(basename $@))).unweighted_unifrac_distance_matrix.qza --o-weighted-unifrac-distance-matrix $(addprefix /,$(basename $(basename $@))).weighted_unifrac_distance_matrix.qza --o-jaccard-distance-matrix $(addprefix /,$(basename $(basename $@))).jaccard_distance_matrix.qza --o-bray-curtis-distance-matrix $(addprefix /,$(basename $(basename $@))).bray_curtis_distance_matrix.qza --o-unweighted-unifrac-pcoa-results $(addprefix /,$(basename $(basename $@))).unweighted_unifrac_pcoa_results.qza --o-weighted-unifrac-pcoa-results $(addprefix /,$(basename $(basename $@))).weighted_unifrac_pcoa_results.qza --o-unweighted-unifrac-emperor $(addprefix /,$(basename $(basename $@))).unweighted_unifrac_emperor.qzv --o-weighted-unifrac-emperor $(addprefix /,$(basename $(basename $@))).weighted_unifrac_emperor.qzv --o-jaccard-emperor $(addprefix /,$(basename $(basename $@))).jaccard_emperor.qzv --o-bray-curtis-emperor $(addprefix /,$(basename $(basename $@))).bray_curtis_emperor.qzv --o-jaccard-pcoa-results $(addprefix /,$(basename $(basename $@))).jaccard_pcoa_results.qza --o-bray-curtis-pcoa-results $(addprefix /,$(basename $(basename $@))).bray_curtis_pcoa_results.qza --p-n-jobs-or-threads $(CPUS) 1> $@.stdout 2> $@.stderr

step29-depth: Output/Step29/everything.DADA2.table.qzv Output/Step29/everything.DADA2.seq.qzv Output/Step29/everything.Deblur.table.qzv Output/Step29/everything.Deblur.seq.qzv
step29: Output/Step29/DADA2/everything.faith_pd_vector.qza Output/Step29/Deblur/everything.faith_pd_vector.qza
.PHONY: step29-depth step29

# Step 30 (Alpha-diversity)
Output/Step30 Output/Step30/DADA2 Output/Step30/Deblur:
	mkdir -p $@

Output/Step30/DADA2/%.tsv: Output/Step29/DADA2/%_vector.qza | Output/Step30/DADA2 build.log
	$(DOCKER) qiime tools export --input-path $(addprefix /,$(word 1,$^)) --output-path $(addprefix /,$(@D)) 1> $@.stdout 2> $@.stderr
	$(DOCKER) mv -f $(addprefix /,$(@D))/alpha-diversity.tsv $(addprefix /,$@)

Output/Step30/Deblur/%.tsv: Output/Step29/Deblur/%_vector.qza | Output/Step30/Deblur build.log
	$(DOCKER) qiime tools export --input-path $(addprefix /,$(word 1,$^)) --output-path $(addprefix /,$(@D)) 1> $@.stdout 2> $@.stderr
	$(DOCKER) mv -f $(addprefix /,$(@D))/alpha-diversity.tsv $(addprefix /,$@)

Output/Step30/everything.DADA2.tar: Python/step30.py Output/Step30/DADA2/everything.evenness.tsv Output/Step30/DADA2/everything.faith_pd.tsv Output/Step30/DADA2/everything.observed_features.tsv Output/Step30/DADA2/everything.shannon.tsv Output/Metadata/everything.tsv | build.log
	$(DOCKER) python3 $(addprefix /,$^ $@) --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/Step30/First.DADA2.tar: Python/step30.py Output/Step30/DADA2/everything.evenness.tsv Output/Step30/DADA2/everything.faith_pd.tsv Output/Step30/DADA2/everything.observed_features.tsv Output/Step30/DADA2/everything.shannon.tsv Output/Metadata/everything.tsv | build.log
	$(DOCKER) python3 $(addprefix /,$^ $@) --cpus $(CPUS) --first 1> $@.stdout 2> $@.stderr

Output/Step30/Second_Third.DADA2.tar: Python/step30.py Output/Step30/DADA2/everything.evenness.tsv Output/Step30/DADA2/everything.faith_pd.tsv Output/Step30/DADA2/everything.observed_features.tsv Output/Step30/DADA2/everything.shannon.tsv Output/Metadata/everything.tsv | build.log
	$(DOCKER) python3 $(addprefix /,$^ $@) --cpus $(CPUS) --second 1> $@.stdout 2> $@.stderr

Output/Step30/everything.Deblur.tar: Python/step30.py Output/Step30/Deblur/everything.evenness.tsv Output/Step30/Deblur/everything.faith_pd.tsv Output/Step30/Deblur/everything.observed_features.tsv Output/Step30/Deblur/everything.shannon.tsv Output/Metadata/everything.tsv | build.log
	$(DOCKER) python3 $(addprefix /,$^ $@) --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/Step30/First.Deblur.tar: Python/step30.py Output/Step30/Deblur/everything.evenness.tsv Output/Step30/Deblur/everything.faith_pd.tsv Output/Step30/Deblur/everything.observed_features.tsv Output/Step30/Deblur/everything.shannon.tsv Output/Metadata/everything.tsv | build.log
	$(DOCKER) python3 $(addprefix /,$^ $@) --cpus $(CPUS) --first 1> $@.stdout 2> $@.stderr

Output/Step30/Second_Third.Deblur.tar: Python/step30.py Output/Step30/Deblur/everything.evenness.tsv Output/Step30/Deblur/everything.faith_pd.tsv Output/Step30/Deblur/everything.observed_features.tsv Output/Step30/Deblur/everything.shannon.tsv Output/Metadata/everything.tsv | build.log
	$(DOCKER) python3 $(addprefix /,$^ $@) --cpus $(CPUS) --second 1> $@.stdout 2> $@.stderr

Output/Step30/Macrogen.tar: Python/step30-1.py Output/TSV/Macrogen.tsv Output/Metadata/everything.tsv | Output/Step30 build.log
	$(DOCKER) python3 $(addprefix /,$^ $@) --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/Step30/First.Macrogen.tar: Python/step30-1.py Output/TSV/Macrogen.tsv Output/Metadata/everything.tsv | Output/Step30 build.log
	$(DOCKER) python3 $(addprefix /,$^ $@) --cpus $(CPUS) --first 1> $@.stdout 2> $@.stderr

Output/Step30/Second_Third.Macrogen.tar: Python/step30-1.py Output/TSV/Macrogen.tsv Output/Metadata/everything.tsv | Output/Step30 build.log
	$(DOCKER) python3 $(addprefix /,$^ $@) --cpus $(CPUS) --second 1> $@.stdout 2> $@.stderr

step30-DADA2: Output/Step30/DADA2/everything.evenness.tsv Output/Step30/DADA2/everything.faith_pd.tsv Output/Step30/DADA2/everything.observed_features.tsv Output/Step30/DADA2/everything.shannon.tsv Output/Step30/everything.DADA2.tar Output/Step30/First.DADA2.tar Output/Step30/Second_Third.DADA2.tar
step30-Deblur: Output/Step30/Deblur/everything.evenness.tsv Output/Step30/Deblur/everything.faith_pd.tsv Output/Step30/Deblur/everything.observed_features.tsv Output/Step30/Deblur/everything.shannon.tsv Output/Step30/everything.Deblur.tar Output/Step30/First.Deblur.tar Output/Step30/Second_Third.Deblur.tar
step30-Macrogen: Output/Step30/Macrogen.tar Output/Step30/First.Macrogen.tar Output/Step30/Second_Third.Macrogen.tar
step30: step30-DADA2 step30-Deblur step30-Macrogen
.PHONY: step30-DADA2 step30-Deblur step30-Macrogen step30

# Step 31 (Violin plots from ANCOM)
Output/ANCOM/everything.DADA2.gg.%.tar: Python/step31.py Output/ANCOM/everything.DADA2.gg.%/ancom.tsv Output/TSV/everything.DADA2.gg.tsv Output/Metadata/everything.tsv | build.log
	$(DOCKER) python3 $(addprefix /,$^ $@) 1> $@.stdout 2> $@.stderr

step31: Output/ANCOM/everything.DADA2.gg.Mouth.pseudocount.Detail_PTB.tar Output/ANCOM/everything.DADA2.gg.Cervix.pseudocount.PTB.tar Output/ANCOM/everything.DADA2.gg.Vagina.pseudocount.PTB.tar Output/ANCOM/everything.DADA2.gg.Baby.pseudocount.PTB.tar Output/ANCOM/everything.DADA2.gg.Mouth.pseudocount.Detail_PTB.tar Output/ANCOM/everything.DADA2.gg.Cervix.pseudocount.Detail_PTB.tar Output/ANCOM/everything.DADA2.gg.Vagina.pseudocount.Detail_PTB.tar Output/ANCOM/everything.DADA2.gg.Baby.pseudocount.Detail_PTB.tar
.PHONY: step31

# Step 32 (Beta-diversity)
Output/Step32/DADA2 Output/Step32/Deblur Output/Step32/Macrogen:
	mkdir -p $@

Output/Step32/DADA2/%.tsv: Output/Step29/DADA2/%_distance_matrix.qza | Output/Step32/DADA2 build.log
	$(DOCKER) qiime tools export --input-path $(addprefix /,$(word 1,$^)) --output-path $(addprefix /,$(@D)) 1> $@.stdout 2> $@.stderr
	$(DOCKER) mv -f $(addprefix /,$(@D))/distance-matrix.tsv $(addprefix /,$@)

Output/Step32/Deblur/%.tsv: Output/Step29/Deblur/%_distance_matrix.qza | Output/Step32/Deblur build.log
	$(DOCKER) qiime tools export --input-path $(addprefix /,$(word 1,$^)) --output-path $(addprefix /,$(@D)) 1> $@.stdout 2> $@.stderr
	$(DOCKER) mv -f $(addprefix /,$(@D))/distance-matrix.tsv $(addprefix /,$@)

Output/Step32/DADA2/%.tar: Python/step32.py Output/Step32/DADA2/%.tsv Output/Metadata/everything.tsv | build.log
	$(DOCKER) python3 $(addprefix /,$^ $@) --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/Step32/DADA2/First.%.tar: Python/step32.py Output/Step32/DADA2/everything.%.tsv Output/Metadata/everything.tsv | build.log
	$(DOCKER) python3 $(addprefix /,$^ $@) --cpus $(CPUS) --first 1> $@.stdout 2> $@.stderr

Output/Step32/DADA2/Second_Third.%.tar: Python/step32.py Output/Step32/DADA2/everything.%.tsv Output/Metadata/everything.tsv | build.log
	$(DOCKER) python3 $(addprefix /,$^ $@) --cpus $(CPUS) --second 1> $@.stdout 2> $@.stderr

Output/Step32/Deblur/%.tar: Python/step32.py Output/Step32/Deblur/%.tsv Output/Metadata/everything.tsv | build.log
	$(DOCKER) python3 $(addprefix /,$^ $@) --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/Step32/Deblur/First.%.tar: Python/step32.py Output/Step32/Deblur/everything.%.tsv Output/Metadata/everything.tsv | build.log
	$(DOCKER) python3 $(addprefix /,$^ $@) --cpus $(CPUS) --first 1> $@.stdout 2> $@.stderr

Output/Step32/Deblur/Second_Third.%.tar: Python/step32.py Output/Step32/Deblur/everything.%.tsv Output/Metadata/everything.tsv | build.log
	$(DOCKER) python3 $(addprefix /,$^ $@) --cpus $(CPUS) --second 1> $@.stdout 2> $@.stderr

Output/Step32/Macrogen/everything.%.tar: Python/step32-1.py Output/TSV/Macrogen.tsv Output/Metadata/everything.tsv | Output/Step32/Macrogen build.log
	$(DOCKER) python3 $(addprefix /,$^ $@) --cpus $(CPUS) --beta "$*" 1> $@.stdout 2> $@.stderr

Output/Step32/Macrogen/First.%.tar: Python/step32-1.py Output/TSV/Macrogen.tsv Output/Metadata/everything.tsv | Output/Step32/Macrogen build.log
	$(DOCKER) python3 $(addprefix /,$^ $@) --cpus $(CPUS) --beta "$*" --first 1> $@.stdout 2> $@.stderr

Output/Step32/Macrogen/Second_Third.%.tar: Python/step32-1.py Output/TSV/Macrogen.tsv Output/Metadata/everything.tsv | Output/Step32/Macrogen build.log
	$(DOCKER) python3 $(addprefix /,$^ $@) --cpus $(CPUS) --beta "$*" --second 1> $@.stdout 2> $@.stderr

step32-DADA2: Output/Step32/DADA2/everything.bray_curtis.tsv Output/Step32/DADA2/everything.bray_curtis.tar Output/Step32/DADA2/First.bray_curtis.tar Output/Step32/DADA2/Second_Third.bray_curtis.tar Output/Step32/DADA2/everything.jaccard.tsv Output/Step32/DADA2/everything.jaccard.tar Output/Step32/DADA2/First.jaccard.tar Output/Step32/DADA2/Second_Third.jaccard.tar Output/Step32/DADA2/everything.unweighted_unifrac.tsv Output/Step32/DADA2/everything.unweighted_unifrac.tar Output/Step32/DADA2/First.unweighted_unifrac.tar Output/Step32/DADA2/Second_Third.unweighted_unifrac.tar Output/Step32/DADA2/everything.weighted_unifrac.tsv Output/Step32/DADA2/everything.weighted_unifrac.tar Output/Step32/DADA2/First.weighted_unifrac.tar Output/Step32/DADA2/Second_Third.weighted_unifrac.tar
step32-Deblur: Output/Step32/Deblur/everything.bray_curtis.tsv Output/Step32/Deblur/everything.bray_curtis.tar Output/Step32/Deblur/First.bray_curtis.tar Output/Step32/Deblur/Second_Third.bray_curtis.tar Output/Step32/Deblur/everything.jaccard.tsv Output/Step32/Deblur/everything.jaccard.tar Output/Step32/Deblur/First.jaccard.tar Output/Step32/Deblur/Second_Third.jaccard.tar Output/Step32/Deblur/everything.unweighted_unifrac.tsv Output/Step32/Deblur/everything.unweighted_unifrac.tar Output/Step32/Deblur/First.unweighted_unifrac.tar Output/Step32/Deblur/Second_Third.unweighted_unifrac.tar Output/Step32/Deblur/everything.weighted_unifrac.tsv Output/Step32/Deblur/everything.weighted_unifrac.tar Output/Step32/Deblur/First.weighted_unifrac.tar Output/Step32/Deblur/Second_Third.weighted_unifrac.tar
step32-Macrogen: Output/Step32/Macrogen/everything.unweighted_unifrac.tar Output/Step32/Macrogen/First.unweighted_unifrac.tar Output/Step32/Macrogen/Second_Third.unweighted_unifrac.tar Output/Step32/Macrogen/everything.weighted_unifrac.tar Output/Step32/Macrogen/First.weighted_unifrac.tar Output/Step32/Macrogen/Second_Third.weighted_unifrac.tar
step32: step32-DADA2 step32-Deblur step32-Macrogen
.PHONY: step32-DADA2 step32-Deblur step32-Macrogen step32

# Step 33 (LEfSe analysis)
Output/Step33 Output/Step33/Default Output/Step33/Cladogram:
	mkdir -p $@
	chmod 777 $@

Output/Step33/%.DADA2.gg.Cervix.tsv Output/Step33/%.DADA2.gg.Mouth.tsv Output/Step33/%.DADA2.gg.Neonate-3day.tsv Output/Step33/%.DADA2.gg.Neonate-5day.tsv Output/Step33/%.DADA2.gg.Vagina.tsv: Python/step33.py Output/TSV/%.DADA2.gg.tsv Output/Metadata/%.tsv | Output/Step33 build.log
	$(DOCKER) python3 $(addprefix /,$^ $(basename $(basename $@))) 1> $@.stdout 2> $@.stderr

Output/Step33/%.DADA2.silva.Cervix.tsv Output/Step33/%.DADA2.silva.Mouth.tsv Output/Step33/%.DADA2.silva.Neonate-3day.tsv Output/Step33/%.DADA2.silva.Neonate-5day.tsv Output/Step33/%.DADA2.silva.Vagina.tsv: Python/step33.py Output/TSV/%.DADA2.silva.tsv Output/Metadata/%.tsv | Output/Step33 build.log
	$(DOCKER) python3 $(addprefix /,$^ $(basename $(basename $@))) 1> $@.stdout 2> $@.stderr

Output/Step33/%.DADA2.homd.Cervix.tsv Output/Step33/%.DADA2.homd.Mouth.tsv Output/Step33/%.DADA2.homd.Neonate-3day.tsv Output/Step33/%.DADA2.homd.Neonate-5day.tsv Output/Step33/%.DADA2.homd.Vagina.tsv: Python/step33.py Output/TSV/%.DADA2.homd.tsv Output/Metadata/%.tsv | Output/Step33 build.log
	$(DOCKER) python3 $(addprefix /,$^ $(basename $(basename $@))) 1> $@.stdout 2> $@.stderr

Output/Step33/%.Deblur.gg.Cervix.tsv Output/Step33/%.Deblur.gg.Mouth.tsv Output/Step33/%.Deblur.gg.Neonate-3day.tsv Output/Step33/%.Deblur.gg.Neonate-5day.tsv Output/Step33/%.Deblur.gg.Vagina.tsv: Python/step33.py Output/TSV/%.Deblur.gg.tsv Output/Metadata/%.tsv | Output/Step33 build.log
	$(DOCKER) python3 $(addprefix /,$^ $(basename $(basename $@))) 1> $@.stdout 2> $@.stderr

Output/Step33/%.Deblur.silva.Cervix.tsv Output/Step33/%.Deblur.silva.Mouth.tsv Output/Step33/%.Deblur.silva.Neonate-3day.tsv Output/Step33/%.Deblur.silva.Neonate-5day.tsv Output/Step33/%.Deblur.silva.Vagina.tsv: Python/step33.py Output/TSV/%.Deblur.silva.tsv Output/Metadata/%.tsv | Output/Step33 build.log
	$(DOCKER) python3 $(addprefix /,$^ $(basename $(basename $@))) 1> $@.stdout 2> $@.stderr

Output/Step33/%.Deblur.homd.Cervix.tsv Output/Step33/%.Deblur.homd.Mouth.tsv Output/Step33/%.Deblur.homd.Neonate-3day.tsv Output/Step33/%.Deblur.homd.Neonate-5day.tsv Output/Step33/%.Deblur.homd.Vagina.tsv: Python/step33.py Output/TSV/%.Deblur.homd.tsv Output/Metadata/%.tsv | Output/Step33 build.log
	$(DOCKER) python3 $(addprefix /,$^ $(basename $(basename $@))) 1> $@.stdout 2> $@.stderr

Output/Step33/%.in: Output/Step33/%.tsv | build.log
	docker run --rm $(VOLUME_OPTS) $(RUN_OPTS) biobakery/lefse format_input.py $(addprefix /,$^ $@) -c 2 -u 1 1> $@.stdout 2> $@.stderr

Output/Step33/%.res: Output/Step33/%.in | build.log
	docker run --rm $(VOLUME_OPTS) $(RUN_OPTS) biobakery/lefse run_lefse.py $(addprefix /,$^ $@) --verbose 1 1> $@.stdout 2> $@.stderr

Output/Step33/Default/%.pdf: Output/Step33/%.res | build.log Output/Step33/Default
	docker run --rm $(VOLUME_OPTS) $(RUN_OPTS) biobakery/lefse plot_res.py $(addprefix /,$^ $@) --format pdf 1> $@.stdout 2> $@.stderr

Output/Step33/Cladogram/%.res: Output/Step33/%.res | build.log Output/Step33/Cladogram
	grep --perl-regexp "True|False" $^ > $@ || cp -f $^ $@

Output/Step33/Cladogram/%.pdf: Output/Step33/Cladogram/%.res | build.log
	docker run --rm $(VOLUME_OPTS) $(RUN_OPTS) biobakery/lefse plot_cladogram.py $(addprefix /,$^ $@) --format pdf 1> $@.stdout 2> $@.stderr

step33-tsv: Output/Step33/everything.DADA2.gg.Cervix.tsv Output/Step33/everything.DADA2.silva.Cervix.tsv Output/Step33/everything.DADA2.homd.Cervix.tsv Output/Step33/everything.Deblur.gg.Cervix.tsv Output/Step33/everything.Deblur.silva.Cervix.tsv Output/Step33/everything.Deblur.homd.Cervix.tsv
step33-in: $(addprefix Output/Step33/everything.DADA2.gg.,$(addsuffix .in,Cervix Mouth Neonate-3day Neonate-5day Vagina)) $(addprefix Output/Step33/everything.DADA2.silva.,$(addsuffix .in,Cervix Mouth Neonate-3day Neonate-5day Vagina)) $(addprefix Output/Step33/everything.DADA2.homd.,$(addsuffix .in,Cervix Mouth Neonate-3day Neonate-5day Vagina)) $(addprefix Output/Step33/everything.Deblur.gg.,$(addsuffix .in,Cervix Mouth Neonate-3day Neonate-5day Vagina)) $(addprefix Output/Step33/everything.Deblur.silva.,$(addsuffix .in,Cervix Mouth Neonate-3day Neonate-5day Vagina)) $(addprefix Output/Step33/everything.Deblur.homd.,$(addsuffix .in,Cervix Mouth Neonate-3day Neonate-5day Vagina))
step33-res: $(addprefix Output/Step33/everything.DADA2.gg.,$(addsuffix .res,Cervix Mouth Neonate-3day Neonate-5day Vagina)) $(addprefix Output/Step33/everything.DADA2.silva.,$(addsuffix .res,Cervix Mouth Neonate-3day Neonate-5day Vagina)) $(addprefix Output/Step33/everything.DADA2.homd.,$(addsuffix .res,Cervix Mouth Neonate-3day Neonate-5day Vagina)) $(addprefix Output/Step33/everything.Deblur.gg.,$(addsuffix .res,Cervix Mouth Neonate-3day Neonate-5day Vagina)) $(addprefix Output/Step33/everything.Deblur.silva.,$(addsuffix .res,Cervix Mouth Neonate-3day Neonate-5day Vagina)) $(addprefix Output/Step33/everything.Deblur.homd.,$(addsuffix .res,Cervix Mouth Neonate-3day Neonate-5day Vagina))
step33-default: $(addprefix Output/Step33/Default/everything.DADA2.gg.,$(addsuffix .pdf,Cervix Mouth Neonate-3day Neonate-5day Vagina)) $(addprefix Output/Step33/Default/everything.DADA2.silva.,$(addsuffix .pdf,Cervix Mouth Neonate-3day Neonate-5day Vagina)) $(addprefix Output/Step33/Default/everything.DADA2.homd.,$(addsuffix .pdf,Cervix Mouth Neonate-3day Neonate-5day Vagina)) $(addprefix Output/Step33/Default/everything.Deblur.gg.,$(addsuffix .pdf,Cervix Mouth Neonate-3day Neonate-5day Vagina)) $(addprefix Output/Step33/Default/everything.Deblur.silva.,$(addsuffix .pdf,Cervix Mouth Neonate-3day Neonate-5day Vagina)) $(addprefix Output/Step33/Default/everything.Deblur.homd.,$(addsuffix .pdf,Cervix Mouth Neonate-3day Neonate-5day Vagina))
step33-cladogram: $(addprefix Output/Step33/Cladogram/everything.DADA2.gg.,$(addsuffix .res,Cervix Mouth Neonate-3day Neonate-5day Vagina)) $(addprefix Output/Step33/Cladogram/everything.DADA2.silva.,$(addsuffix .res,Cervix Mouth Neonate-3day Neonate-5day Vagina)) $(addprefix Output/Step33/Cladogram/everything.DADA2.homd.,$(addsuffix .res,Cervix Mouth Neonate-3day Neonate-5day Vagina)) $(addprefix Output/Step33/Cladogram/everything.Deblur.gg.,$(addsuffix .res,Cervix Mouth Neonate-3day Neonate-5day Vagina)) $(addprefix Output/Step33/Cladogram/everything.Deblur.silva.,$(addsuffix .res,Cervix Mouth Neonate-3day Neonate-5day Vagina)) $(addprefix Output/Step33/Cladogram/everything.Deblur.homd.,$(addsuffix .res,Cervix Mouth Neonate-3day Neonate-5day Vagina)) $(addprefix Output/Step33/Cladogram/everything.DADA2.gg.,$(addsuffix .pdf,Cervix Mouth Neonate-3day Neonate-5day Vagina)) $(addprefix Output/Step33/Cladogram/everything.DADA2.silva.,$(addsuffix .pdf,Cervix Mouth Neonate-3day Neonate-5day Vagina)) $(addprefix Output/Step33/Cladogram/everything.DADA2.homd.,$(addsuffix .pdf,Cervix Mouth Neonate-3day Neonate-5day Vagina)) $(addprefix Output/Step33/Cladogram/everything.Deblur.gg.,$(addsuffix .pdf,Cervix Mouth Neonate-3day Neonate-5day Vagina)) $(addprefix Output/Step33/Cladogram/everything.Deblur.silva.,$(addsuffix .pdf,Cervix Mouth Neonate-3day Neonate-5day Vagina)) $(addprefix Output/Step33/Cladogram/everything.Deblur.homd.,$(addsuffix .pdf,Cervix Mouth Neonate-3day Neonate-5day Vagina))
step33: step33-tsv step33-in step33-res step33-default step33-cladogram
.PHONY: step33-tsv step33-in step33-res step33-default step33-cladogram step33

# Step 34 (LEfSe analysis - Obesity)
Output/Step34 Output/Step34/Default Output/Step34/Cladogram:
	mkdir -p $@
	chmod 777 $@

Output/Step34/%.DADA2.gg.Cervix.tsv Output/Step34/%.DADA2.gg.Mouth.tsv Output/Step34/%.DADA2.gg.Neonate-3day.tsv Output/Step34/%.DADA2.gg.Neonate-5day.tsv Output/Step34/%.DADA2.gg.Vagina.tsv: Python/step33.py Output/TSV/%.DADA2.gg.tsv Output/Metadata/%.tsv | Output/Step34 build.log
	$(DOCKER) python3 $(addprefix /,$^ $(basename $(basename $@))) --c "Obesity" 1> $@.stdout 2> $@.stderr

Output/Step34/%.DADA2.silva.Cervix.tsv Output/Step34/%.DADA2.silva.Mouth.tsv Output/Step34/%.DADA2.silva.Neonate-3day.tsv Output/Step34/%.DADA2.silva.Neonate-5day.tsv Output/Step34/%.DADA2.silva.Vagina.tsv: Python/step33.py Output/TSV/%.DADA2.silva.tsv Output/Metadata/%.tsv | Output/Step34 build.log
	$(DOCKER) python3 $(addprefix /,$^ $(basename $(basename $@))) --c "Obesity" 1> $@.stdout 2> $@.stderr

Output/Step34/%.DADA2.homd.Cervix.tsv Output/Step34/%.DADA2.homd.Mouth.tsv Output/Step34/%.DADA2.homd.Neonate-3day.tsv Output/Step34/%.DADA2.homd.Neonate-5day.tsv Output/Step34/%.DADA2.homd.Vagina.tsv: Python/step33.py Output/TSV/%.DADA2.homd.tsv Output/Metadata/%.tsv | Output/Step34 build.log
	$(DOCKER) python3 $(addprefix /,$^ $(basename $(basename $@))) --c "Obesity" 1> $@.stdout 2> $@.stderr

Output/Step34/%.Deblur.gg.Cervix.tsv Output/Step34/%.Deblur.gg.Mouth.tsv Output/Step34/%.Deblur.gg.Neonate-3day.tsv Output/Step34/%.Deblur.gg.Neonate-5day.tsv Output/Step34/%.Deblur.gg.Vagina.tsv: Python/step33.py Output/TSV/%.Deblur.gg.tsv Output/Metadata/%.tsv | Output/Step34 build.log
	$(DOCKER) python3 $(addprefix /,$^ $(basename $(basename $@))) --c "Obesity" 1> $@.stdout 2> $@.stderr

Output/Step34/%.Deblur.silva.Cervix.tsv Output/Step34/%.Deblur.silva.Mouth.tsv Output/Step34/%.Deblur.silva.Neonate-3day.tsv Output/Step34/%.Deblur.silva.Neonate-5day.tsv Output/Step34/%.Deblur.silva.Vagina.tsv: Python/step33.py Output/TSV/%.Deblur.silva.tsv Output/Metadata/%.tsv | Output/Step34 build.log
	$(DOCKER) python3 $(addprefix /,$^ $(basename $(basename $@))) --c "Obesity" 1> $@.stdout 2> $@.stderr

Output/Step34/%.Deblur.homd.Cervix.tsv Output/Step34/%.Deblur.homd.Mouth.tsv Output/Step34/%.Deblur.homd.Neonate-3day.tsv Output/Step34/%.Deblur.homd.Neonate-5day.tsv Output/Step34/%.Deblur.homd.Vagina.tsv: Python/step33.py Output/TSV/%.Deblur.homd.tsv Output/Metadata/%.tsv | Output/Step34 build.log
	$(DOCKER) python3 $(addprefix /,$^ $(basename $(basename $@))) --c "Obesity" 1> $@.stdout 2> $@.stderr

Output/Step34/%.Macrogen.Cervix.tsv Output/Step34/%.Macrogen.Mouth.tsv Output/Step34/%.Macrogen.Neonate-3day.tsv Output/Step34/%.Macrogen.Neonate-5day.tsv Output/Step34/%.Macrogen.Vagina.tsv: Python/step33.py Output/TSV/Macrogen.tsv Output/Metadata/%.tsv | Output/Step34 build.log
	$(DOCKER) python3 $(addprefix /,$^ $(basename $(basename $@))) --c "Obesity" 1> $@.stdout 2> $@.stderr

Output/Step34/%.in: Output/Step34/%.tsv | build.log
	docker run --rm $(VOLUME_OPTS) $(RUN_OPTS) biobakery/lefse format_input.py $(addprefix /,$^ $@) -c 2 -u 1 1> $@.stdout 2> $@.stderr

Output/Step34/%.res: Output/Step34/%.in | build.log
	docker run --rm $(VOLUME_OPTS) $(RUN_OPTS) biobakery/lefse run_lefse.py $(addprefix /,$^ $@) --verbose 1 1> $@.stdout 2> $@.stderr || touch $@

Output/Step34/Default/%.pdf: Output/Step34/%.res | build.log Output/Step34/Default
	docker run --rm $(VOLUME_OPTS) $(RUN_OPTS) biobakery/lefse plot_res.py $(addprefix /,$^ $@) --format pdf 1> $@.stdout 2> $@.stderr || touch $@

Output/Step34/Cladogram/%.res: Output/Step34/%.res | build.log Output/Step34/Cladogram
	grep --perl-regexp "True|False" $^ > $@ || cp -f $^ $@

Output/Step34/Cladogram/%.pdf: Output/Step34/Cladogram/%.res | build.log
	docker run --rm $(VOLUME_OPTS) $(RUN_OPTS) biobakery/lefse plot_cladogram.py $(addprefix /,$^ $@) --format pdf 1> $@.stdout 2> $@.stderr || touch $@

step34-tsv: Output/Step34/everything.DADA2.gg.Cervix.tsv Output/Step34/everything.DADA2.silva.Cervix.tsv Output/Step34/everything.DADA2.homd.Cervix.tsv Output/Step34/everything.Deblur.gg.Cervix.tsv Output/Step34/everything.Deblur.silva.Cervix.tsv Output/Step34/everything.Deblur.homd.Cervix.tsv Output/Step34/everything.Macrogen.Cervix.tsv
step34-in: $(addprefix Output/Step34/everything.DADA2.gg.,$(addsuffix .in,Cervix Mouth Neonate-3day Neonate-5day Vagina)) $(addprefix Output/Step34/everything.DADA2.silva.,$(addsuffix .in,Cervix Mouth Neonate-3day Neonate-5day Vagina)) $(addprefix Output/Step34/everything.DADA2.homd.,$(addsuffix .in,Cervix Mouth Neonate-3day Neonate-5day Vagina)) $(addprefix Output/Step34/everything.Deblur.gg.,$(addsuffix .in,Cervix Mouth Neonate-3day Neonate-5day Vagina)) $(addprefix Output/Step34/everything.Deblur.silva.,$(addsuffix .in,Cervix Mouth Neonate-3day Neonate-5day Vagina)) $(addprefix Output/Step34/everything.Deblur.homd.,$(addsuffix .in,Cervix Mouth Neonate-3day Neonate-5day Vagina)) $(addprefix Output/Step34/everything.Macrogen.,$(addsuffix .in,Cervix Mouth Neonate-3day Neonate-5day Vagina))
step34-res: $(addprefix Output/Step34/everything.DADA2.gg.,$(addsuffix .res,Cervix Mouth Neonate-3day Neonate-5day Vagina)) $(addprefix Output/Step34/everything.DADA2.silva.,$(addsuffix .res,Cervix Mouth Neonate-3day Neonate-5day Vagina)) $(addprefix Output/Step34/everything.DADA2.homd.,$(addsuffix .res,Cervix Mouth Neonate-3day Neonate-5day Vagina)) $(addprefix Output/Step34/everything.Deblur.gg.,$(addsuffix .res,Cervix Mouth Neonate-3day Neonate-5day Vagina)) $(addprefix Output/Step34/everything.Deblur.silva.,$(addsuffix .res,Cervix Mouth Neonate-3day Neonate-5day Vagina)) $(addprefix Output/Step34/everything.Deblur.homd.,$(addsuffix .res,Cervix Mouth Neonate-3day Neonate-5day Vagina)) $(addprefix Output/Step34/everything.Macrogen.,$(addsuffix .res,Cervix Mouth Neonate-3day Neonate-5day Vagina))
step34-default: $(addprefix Output/Step34/Default/everything.DADA2.gg.,$(addsuffix .pdf,Cervix Mouth Neonate-3day Neonate-5day Vagina)) $(addprefix Output/Step34/Default/everything.DADA2.silva.,$(addsuffix .pdf,Cervix Mouth Neonate-3day Neonate-5day Vagina)) $(addprefix Output/Step34/Default/everything.DADA2.homd.,$(addsuffix .pdf,Cervix Mouth Neonate-3day Neonate-5day Vagina)) $(addprefix Output/Step34/Default/everything.Deblur.gg.,$(addsuffix .pdf,Cervix Mouth Neonate-3day Neonate-5day Vagina)) $(addprefix Output/Step34/Default/everything.Deblur.silva.,$(addsuffix .pdf,Cervix Mouth Neonate-3day Neonate-5day Vagina)) $(addprefix Output/Step34/Default/everything.Deblur.homd.,$(addsuffix .pdf,Cervix Mouth Neonate-3day Neonate-5day Vagina)) $(addprefix Output/Step34/Default/everything.Macrogen.,$(addsuffix .pdf,Cervix Mouth Neonate-3day Neonate-5day Vagina))
step34-cladogram: $(addprefix Output/Step34/Cladogram/everything.DADA2.gg.,$(addsuffix .res,Cervix Mouth Neonate-3day Neonate-5day Vagina)) $(addprefix Output/Step34/Cladogram/everything.DADA2.silva.,$(addsuffix .res,Cervix Mouth Neonate-3day Neonate-5day Vagina)) $(addprefix Output/Step34/Cladogram/everything.DADA2.homd.,$(addsuffix .res,Cervix Mouth Neonate-3day Neonate-5day Vagina)) $(addprefix Output/Step34/Cladogram/everything.Deblur.gg.,$(addsuffix .res,Cervix Mouth Neonate-3day Neonate-5day Vagina)) $(addprefix Output/Step34/Cladogram/everything.Deblur.silva.,$(addsuffix .res,Cervix Mouth Neonate-3day Neonate-5day Vagina)) $(addprefix Output/Step34/Cladogram/everything.Deblur.homd.,$(addsuffix .res,Cervix Mouth Neonate-3day Neonate-5day Vagina)) $(addprefix Output/Step34/Cladogram/everything.DADA2.gg.,$(addsuffix .pdf,Cervix Mouth Neonate-3day Neonate-5day Vagina)) $(addprefix Output/Step34/Cladogram/everything.DADA2.silva.,$(addsuffix .pdf,Cervix Mouth Neonate-3day Neonate-5day Vagina)) $(addprefix Output/Step34/Cladogram/everything.DADA2.homd.,$(addsuffix .pdf,Cervix Mouth Neonate-3day Neonate-5day Vagina)) $(addprefix Output/Step34/Cladogram/everything.Deblur.gg.,$(addsuffix .pdf,Cervix Mouth Neonate-3day Neonate-5day Vagina)) $(addprefix Output/Step34/Cladogram/everything.Deblur.silva.,$(addsuffix .pdf,Cervix Mouth Neonate-3day Neonate-5day Vagina)) $(addprefix Output/Step34/Cladogram/everything.Deblur.homd.,$(addsuffix .pdf,Cervix Mouth Neonate-3day Neonate-5day Vagina)) $(addprefix Output/Step34/Cladogram/everything.Macrogen.,$(addsuffix .pdf,Cervix Mouth Neonate-3day Neonate-5day Vagina))
step34: step34-tsv step34-in step34-res step34-default step34-cladogram
.PHONY: step34-tsv step34-in step34-res step34-default step34-cladogram step34

# Step 35 (Cluster map plot within bacteria)
Output/Step35:
	mkdir -p $@

Output/Step35/%.pearson.png: Python/step35.py Output/TSV/%.tar.gz | build.log Output/Step35
	$(DOCKER) python3 $(addprefix /,$^ $@) --cpus $(CPUS) --pearson 1> $@.stdout 2> $@.stderr

Output/Step35/%.spearman.png: Python/step35.py Output/TSV/%.tar.gz | build.log Output/Step35
	$(DOCKER) python3 $(addprefix /,$^ $@) --cpus $(CPUS) --spearman 1> $@.stdout 2> $@.stderr

Output/Step35/%.kendall.png: Python/step35.py Output/TSV/%.tar.gz | build.log Output/Step35
	$(DOCKER) python3 $(addprefix /,$^ $@) --cpus $(CPUS) --kendall 1> $@.stdout 2> $@.stderr

Output/Step35/Macrogen.%.png: Python/step35-1.py Output/TSV/Macrogen.tar.gz | build.log Output/Step35
	$(DOCKER) python3 $(addprefix /,$^ $@) --cpus $(CPUS) --$* 1> $@.stdout 2> $@.stderr

step35-everything-DADA2: Output/Step35/everything.DADA2.gg.pearson.png Output/Step35/everything.DADA2.gg.spearman.png Output/Step35/everything.DADA2.gg.kendall.png Output/Step35/everything.DADA2.silva.pearson.png Output/Step35/everything.DADA2.silva.spearman.png Output/Step35/everything.DADA2.silva.kendall.png Output/Step35/everything.DADA2.homd.pearson.png Output/Step35/everything.DADA2.homd.spearman.png Output/Step35/everything.DADA2.homd.kendall.png
step35-everything-Deblur: Output/Step35/everything.Deblur.gg.pearson.png Output/Step35/everything.Deblur.gg.spearman.png Output/Step35/everything.Deblur.gg.kendall.png Output/Step35/everything.Deblur.silva.pearson.png Output/Step35/everything.Deblur.silva.spearman.png Output/Step35/everything.Deblur.silva.kendall.png Output/Step35/everything.Deblur.homd.pearson.png Output/Step35/everything.Deblur.homd.spearman.png Output/Step35/everything.Deblur.homd.kendall.png
step35-Stool-DADA2: Output/Step35/Stool.DADA2.gg.pearson.png Output/Step35/Stool.DADA2.gg.spearman.png Output/Step35/Stool.DADA2.gg.kendall.png Output/Step35/Stool.DADA2.silva.pearson.png Output/Step35/Stool.DADA2.silva.spearman.png Output/Step35/Stool.DADA2.silva.kendall.png Output/Step35/Stool.DADA2.homd.pearson.png Output/Step35/Stool.DADA2.homd.spearman.png Output/Step35/Stool.DADA2.homd.kendall.png
step35-Stool-Deblur: Output/Step35/Stool.Deblur.gg.pearson.png Output/Step35/Stool.Deblur.gg.spearman.png Output/Step35/Stool.Deblur.gg.kendall.png Output/Step35/Stool.Deblur.silva.pearson.png Output/Step35/Stool.Deblur.silva.spearman.png Output/Step35/Stool.Deblur.silva.kendall.png Output/Step35/Stool.Deblur.homd.pearson.png Output/Step35/Stool.Deblur.homd.spearman.png Output/Step35/Stool.Deblur.homd.kendall.png
step35-Macrogen: Output/Step35/Macrogen.pearson.png Output/Step35/Macrogen.spearman.png Output/Step35/Macrogen.kendall.png
step35: step35-everything-DADA2 step35-everything-Deblur step35-Stool-DADA2 step35-Stool-Deblur step35-Macrogen
.PHONY: step35-everything-DADA2 step35-everything-Deblur step35-Stool-DADA2 step35-Stool-Deblur step35-Macrogen step35

# Step 36 (Scatter plot heatmap within bacteria)
Output/Step36:
	mkdir -p $@

Output/Step36/%.pearson.png: Python/step36.py Output/TSV/%.tar.gz | build.log Output/Step36
	$(DOCKER) python3 $(addprefix /,$^ $@) --cpus $(CPUS) --pearson 1> $@.stdout 2> $@.stderr

Output/Step36/%.spearman.png: Python/step36.py Output/TSV/%.tar.gz | build.log Output/Step36
	$(DOCKER) python3 $(addprefix /,$^ $@) --cpus $(CPUS) --spearman 1> $@.stdout 2> $@.stderr

Output/Step36/%.kendall.png: Python/step36.py Output/TSV/%.tar.gz | build.log Output/Step36
	$(DOCKER) python3 $(addprefix /,$^ $@) --cpus $(CPUS) --kendall 1> $@.stdout 2> $@.stderr

Output/Step36/Macrogen.%.png: Python/step36-1.py Output/TSV/Macrogen.tar.gz | build.log Output/Step36
	$(DOCKER) python3 $(addprefix /,$^ $@) --$* --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

step36-everything-DADA2: Output/Step36/everything.DADA2.gg.pearson.png Output/Step36/everything.DADA2.gg.spearman.png Output/Step36/everything.DADA2.gg.kendall.png Output/Step36/everything.DADA2.silva.pearson.png Output/Step36/everything.DADA2.silva.spearman.png Output/Step36/everything.DADA2.silva.kendall.png Output/Step36/everything.DADA2.homd.pearson.png Output/Step36/everything.DADA2.homd.spearman.png Output/Step36/everything.DADA2.homd.kendall.png
step36-everything-Deblur: Output/Step36/everything.Deblur.gg.pearson.png Output/Step36/everything.Deblur.gg.spearman.png Output/Step36/everything.Deblur.gg.kendall.png Output/Step36/everything.Deblur.silva.pearson.png Output/Step36/everything.Deblur.silva.spearman.png Output/Step36/everything.Deblur.silva.kendall.png Output/Step36/everything.Deblur.homd.pearson.png Output/Step36/everything.Deblur.homd.spearman.png Output/Step36/everything.Deblur.homd.kendall.png
step36-Stool-DADA2: Output/Step36/Stool.DADA2.gg.pearson.png Output/Step36/Stool.DADA2.gg.spearman.png Output/Step36/Stool.DADA2.gg.kendall.png Output/Step36/Stool.DADA2.silva.pearson.png Output/Step36/Stool.DADA2.silva.spearman.png Output/Step36/Stool.DADA2.silva.kendall.png Output/Step36/Stool.DADA2.homd.pearson.png Output/Step36/Stool.DADA2.homd.spearman.png Output/Step36/Stool.DADA2.homd.kendall.png
step36-Stool-Deblur: Output/Step36/Stool.Deblur.gg.pearson.png Output/Step36/Stool.Deblur.gg.spearman.png Output/Step36/Stool.Deblur.gg.kendall.png Output/Step36/Stool.Deblur.silva.pearson.png Output/Step36/Stool.Deblur.silva.spearman.png Output/Step36/Stool.Deblur.silva.kendall.png Output/Step36/Stool.Deblur.homd.pearson.png Output/Step36/Stool.Deblur.homd.spearman.png Output/Step36/Stool.Deblur.homd.kendall.png
step36-Macrogen: Output/Step36/Macrogen.pearson.png Output/Step36/Macrogen.spearman.png Output/Step36/Macrogen.kendall.png
step36: step36-everything-DADA2 step36-everything-Deblur step36-Stool-DADA2 step36-Stool-Deblur step36-Macrogen
.PHONY: step36-everything-DADA2 step36-everything-Deblur step36-Stool-DADA2 step36-Stool-Deblur step36-Macrogen step36

# Step37 (Correlation within sites)
Output/Step37:
	mkdir -p $@

Output/Step37/Macrogen.%.pdf: Python/step37-1.py Output/TSV/Macrogen.tar.gz | build.log Output/Step37
	$(DOCKER) python3 $(addprefix /,$^ $@) --$* --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

step37-Macrogen: Output/Step37/Macrogen.pearson.pdf Output/Step37/Macrogen.spearman.pdf Output/Step37/Macrogen.kendall.pdf
step37: step37-Macrogen
.PHONY: step37-Macrogen step37

# Step 38 (Violin plots within sites)
Output/Step38:
	mkdir -p $@

Output/Step37/Macrogen.%.pdf: Python/step37-1.py Output/TSV/Macrogen.tar.gz | build.log Output/Step37
	$(DOCKER) python3 $(addprefix /,$^ $@) --$* --cpus $(CPUS)
