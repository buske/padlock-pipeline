trimmomatic_version := 0.32
trimmomatic := Trimmomatic-$(trimmomatic_version)
trimmomatic_bin := lib/$(trimmomatic)/trimmomatic-$(trimmomatic_version).jar

bismark_version := v0.12.3
bismark := bismark_$(bismark_version)
bismark_bin := lib/$(bismark)/bismark

bowtie_version := 1.1.1
bowtie := bowtie-$(bowtie_version)-linux-x86_64
bowtie_bin := lib/bowtie-$(bowtie_version)/bowtie

install: bismark trimmomatic bowtie

bismark: $(bismark_bin)

trimmomatic: $(trimmomatic_bin)

bowtie: $(bowtie_bin)


lib/$(trimmomatic).zip:
	curl -L -o $@ http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/$(trimmomatic).zip

$(trimmomatic_bin): lib/$(trimmomatic).zip
	cd lib; unzip $(<F)
	touch $@

lib/$(bismark).tar.gz:
	curl -L -o $@ http://www.bioinformatics.babraham.ac.uk/projects/bismark/$(bismark).tar.gz

$(bismark_bin): lib/$(bismark).tar.gz
	cd lib; tar -xzf $(<F)
	touch $@


lib/$(bowtie).zip:
	curl -L -o $@ http://downloads.sourceforge.net/project/bowtie-bio/bowtie/1.1.1/$(bowtie).zip

$(bowtie_bin): lib/$(bowtie).zip
	cd lib; unzip $(<F)
	touch $@
