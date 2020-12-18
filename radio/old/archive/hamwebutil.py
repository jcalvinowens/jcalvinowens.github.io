#!/usr/bin/env python3
#
# Copyright (C) 2017 Calvin Owens <jcalvinowens@gmail.com>
#
# This program is free software; you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation; version 2 only.
#
# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
# details.

import argparse
import copy
import datetime
import itertools
import json
import math
import sys
import time

from collections import defaultdict
import matplotlib.pyplot as plt

from adifutil import QSO, QSOGroup

class WobsiteGroup(QSOGroup):
	FIGURESIZE = (10,4)
	BANDCOLORS = [
		(160, 'k'),
		(80, 'k'),
		(40, 'r'),
		(30, 'y'),
		(20, 'g'),
		(17, 'b'),
		(15, 'c'),
		(12, 'm'),
		(10, 'm'),
		(6, 'm'),
	]

	def save_total_pattern(self, filename, rmax=12450):
		def mkt(d):
			hr = (d.dateon.hour + 16) % 24
			s = hr * 3600 + d.dateon.minute * 60 + d.dateon.second
			return s / 86400 * (math.pi * 2)

		fig = plt.figure(tight_layout=True, figsize=(10,10))
		p = plt.subplot(111, projection="polar")
		for band, color in self.BANDCOLORS:
			p.plot(
				[x.azrads for x in self.qsos if x.band == band],
				[x.dist for x in self.qsos if x.band == band], color + "."
			)
		p.set_rticks([3500, 6000, 8000, 10000, rmax] if rmax > 3500 else
			([1000, 2000, 3000, rmax] if rmax > 1000 else [250, 500,
			750, rmax]))
		p.set_rmax(rmax)
		p.set_rmin(0)
		p.set_theta_zero_location("N")
		p.set_theta_direction(-1)
		fig.savefig(filename, format=filename[-3:], bbox_inches='tight')

	def save_band_patterns(self, filename, rmax=12450):
		i = 1

		fig = plt.figure(tight_layout=True, figsize=self.FIGURESIZE)
		for band, color in self.BANDCOLORS:
			p = plt.subplot(2, 5, i, projection="polar")
			p.set_title("%dm Band" % (band))
			p.plot([x.azrads for x in self.qsos if x.band == band],
			       [x.dist for x in self.qsos if x.band == band], color + ".")
			p.set_theta_zero_location("N")
			p.set_theta_direction(-1)
			p.set_rticks([])
			p.set_xticklabels([])
			p.set_rmax(rmax)
			p.set_rmin(0)
			i += 1

		fig.savefig(filename, format=filename[-3:])

	def save_sigstrengthrcvd_plots(self, filename):
		i = 1

		fig = plt.figure(tight_layout=True, figsize=self.FIGURESIZE)
		for band, color in self.BANDCOLORS:
			p = plt.subplot(2, 5, i)
			p.set_title("%dm Band" % (band))
			p.set_xlabel("Distance (mi)")
			p.set_ylabel("RST Received (db)")
			p.axis([100, 12450, -25, 15])
			p.semilogx([x.dist for x in self.qsos if x.band == band],
				   [x.rst_rcvd for x in self.qsos if x.band == band], color + ".")
			i += 1

		fig.savefig(filename, format=filename[-3:])

	def save_sigstrengthrel_plots(self, filename):
		i = 1

		fig = plt.figure(tight_layout=True, figsize=self.FIGURESIZE)
		for band, color in self.BANDCOLORS:
			p = plt.subplot(2, 5, i)
			p.set_title("%dm Band" % (band))
			p.set_xlabel("RST Sent (db)")
			p.set_ylabel("RST Received (db)")
			p.axis([-25, 15, -25, 15])
			p.plot([x.rst_sent for x in self.qsos if x.band == band],
			       [x.rst_rcvd for x in self.qsos if x.band == band], color + ".")
			i += 1

		fig.savefig(filename, format=filename[-3:])

	def save_clock_dist_plot(self, filename):
		def mkt(d):
			hr = (d.dateon.hour + 16) % 24
			s = hr * 3600 + d.dateon.minute * 60 + d.dateon.second
			return s / 86400 * (math.pi * 2)

		fig = plt.figure(tight_layout=True, figsize=(10,10))

		p = plt.subplot(121, projection="polar")
		for band, color in self.BANDCOLORS:
			p.plot(
				[mkt(x) for x in self.qsos if x.band == band],
				[x.dist for x in self.qsos if x.band == band], color + "."
			)
		p.set_rticks([1000, 2000, 4000, 8000])
		p.set_rmax(12450)
		p.set_rmin(0)
		p.set_xticklabels(['00:00', '03:00', '06:00', '09:00', '12:00', '15:00', '18:00', '21:00'])
		p.set_theta_zero_location("N")
		p.set_theta_direction(-1)

		p = plt.subplot(122, projection="polar")
		for band, color in self.BANDCOLORS:
			p.plot(
				[mkt(x) for x in self.qsos if x.band == band],
				[x.dist for x in self.qsos if x.band == band], color + "."
			)
		p.set_rticks([500, 1000, 2000])
		p.set_rmax(3500)
		p.set_rmin(0)
		p.set_xticklabels(['00:00', '03:00', '06:00', '09:00', '12:00', '15:00', '18:00', '21:00'])
		p.set_theta_zero_location("N")
		p.set_theta_direction(-1)

		fig.savefig(filename, format=filename[-3:], bbox_inches='tight')

	def save_all_charts(self, d):
		self.save_total_pattern(filename=f"{d}/all-pattern.svg")
		self.save_total_pattern(filename=f"{d}/all-pattern-close.svg", rmax=3500)
		self.save_total_pattern(filename=f"{d}/all-pattern-closer.svg", rmax=1000)

	def html_table(self, fd, nr=20, srt=lambda k: k.dateon,
		       fil=lambda x: True, caption=None):
		fd.write("<table class=\"report\">")

		if caption:
			fd.write(f"<caption class=\"reportcaption\">{caption}</caption>")

		hdrs = [
			"Callsign",
			"Grid",
			"Date UTC",
			"Hdg",
			"Dist",
			"Band",
			"RST TX/RX",
			"DXCC",
		]

		# Column headers
		fd.write("<tr>")
		for header in hdrs:
			fd.write(f"<th>{header}</th>")
		fd.write("</tr>")

		# Data
		for qso in sorted([x for x in self.qsos if fil(x)],
				  key=srt, reverse=True)[:nr]:
			country = qso.country
			if qso.country == "UNITED STATES OF AMERICA":
				country = f"{qso.state}, USA"

			if country:
				country = country.replace("REPUBLIC OF", "R.O.")

			if not country:
				country = "<i>Unconfirmed</i>"

			cols = [
				qso.call,
				f"{qso.gridsquare[:4]}",
				qso.dateon.isoformat(" ")[:-15],
				qso.cardinal,
				f"{qso.dist}mi",
				f"{qso.band}m",
				f"{qso.rst_sent:+03d}/{qso.rst_rcvd:+03d}" if qso.rst_sent is not None else "N/A",
				country,
			]

			fd.write("<tr>")
			for val in cols:
				fd.write(f"<td>{val}</td>")
			fd.write("</tr>")

		fd.write("</table>")

	def get_states_list(self):
		states = set(qso.state for qso in self.qsos if qso.country == "UNITED STATES OF AMERICA")
		states.update(qso.state for qso in self.qsos if qso.country == "ALASKA")
		states.update(qso.state for qso in self.qsos if qso.country == "HAWAII")
		ret = ""

		for i, state in enumerate(sorted(states)):
			ret += state + "(%02d) " % len([x for x in self.qsos if x.state == state])
			if (i + 1) % 5 == 0:
				ret += "\n"
		ret += "\n"
		return ret

	def get_dxcc_list(self):
		def strfn(x):
			l = [q for q in self.qsos if q.country == x.country]
			if len(l) > 1:
				return f"({len(l)})"

			return ""

		x = set(f"  {x.country} {strfn(x)}"
			for x in self.qsos if x.country)

		return "\n".join(sorted(x))


	def save_stats(self, d):
		with open(f"{d}/stats.html", "w") as fd:
			fd.write("<!doctype html>")
			fd.write(f"""<pre style="column-count: 2">
Generated {datetime.datetime.now().isoformat(" ")}

Totals:

  {len(self.qsos)} QSOs, {len(self.qsls)} Confirmed ({len(self.qsls) / len(self.qsos) * 100:.1f}%)

Confirmed DXCC Countries ({len(set(x.country for x in self.qsos if x.country))}):

{self.get_dxcc_list()}

Confirmed US States:

{self.get_states_list()}

</pre>""")

		with open(f"{d}/tables.html", "w") as fd:
			fd.write("<!doctype html>\n\n")

			self.html_table(fd, nr=25, caption="Most Recent QSOs")
			fd.write("<hr class=\"report\">")
			self.html_table(fd, nr=40, srt=lambda k: k.dist,
					fil=lambda q: q.qsl_rcvd,
					caption="Most Distant Confirmed QSLs")

		with open(f"{d}/all-contacts.html", "w") as fd:
			fd.write("<!doctype html>\n\n")
			self.html_table(fd, nr=99999, caption="All Contacts")


if __name__ == "__main__":
	qsos = WobsiteGroup(QSO.load_qsos_from_files(*sys.argv[1:]))
	home_qsos = WobsiteGroup(qsos.qsos, call="W5CAL", mode="FT8",
				 grid="CM87WL")

	home_qsos.save_all_charts("gen/")
	qsos.save_stats("gen/")

	with open("gen/tag.txt", "w+") as tagfile:
		tagfile.write(str(time.time()))
