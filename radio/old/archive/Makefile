HAMWEBUTIL = ~/git/miscellaneous/radio/hamwebutil.py
ADIFLOGS = ~/.local/share/WSJT-X/wsjtx_log.adi ~/Downloads/lotwreport.adi

pygen = gen/tag.txt
pages = gen/tables.html gen/stats.html
charts = gen/all-pattern.svg gen/all-pattern-close.svg

all: $(pages) $(charts)
$(pages): $(pygen)
$(charts): $(pygen)
$(pygen): $(HAMWEBUTIL) $(ADIFLOGS)

$(pygen):
	$(HAMWEBUTIL) $(ADIFLOGS)

clean:
	rm -f $(pygen) $(pages) $(charts)
