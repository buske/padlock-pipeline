trimmomatic_version := 0.32
trimmomatic := Trimmomatic-$(trimmomatic_version)
trimmomatic_dir := lib/$(trimmomatic)/
bismark_version := v0.12.3
bismark := bismark_$(bismark_version)
bismark_dir := lib/$(bismark)/

install: $(bismark_bin) $(trimmomatic_jar)

lib/$(trimmomatic).zip:
	curl -o $@ http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/$(trimmomatic).zip

$(trimmomatic_dir): lib/$(trimmomatic).zip
	cd lib; unzip $(<F)

lib/$(bismark).tar.gz:
	curl -o $@ http://www.bioinformatics.babraham.ac.uk/projects/bismark/$(bismark).tar.gz

$(bismark_dir): lib/$(bismark).tar.gz
	cd lib; tar -xzf $(<F)
