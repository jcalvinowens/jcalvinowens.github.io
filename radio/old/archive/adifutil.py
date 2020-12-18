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

from collections import defaultdict

class QSO(object):
	FUZZTIME = datetime.timedelta(minutes=30)
	CARDINALS = ["N", "NNE", "NE", "ENE", "E", "ESE", "SE", "SSE", "S",
		     "SSW", "SW", "WSW", "W", "WNW", "NW", "NNW", "N"]
	REQUIRED_FIELDS = ["call", "qso_date", "band"]
	DEFAULT_FIELDS = {
		"qsl_rcvd": False,
	}
	PARSED_FIELDS = {
		"band": lambda v: int("".join([s for s in v if s.isdigit()])),
		"tx_pwr": lambda v: v,
		"qsl_rcvd": lambda v: True if v == "Y" else False,
		"freq": float,
		"rst_sent": int,
		"rst_rcvd": int,
		"my_cq_zone": int,
		"my_itu_zone": int,
		"my_dxcc": int,
		"cqz": int,
		"ituz": int,
		"dxcc": int,
	}

	@staticmethod
	def next_adif_tag(s):
		"""Parse the next ADIF tag in stream @s. Returns the key, the
		value (or None for valueless tags), and the number of characters
		the caller should discard from @s before calling us again."""

		skipped = 0
		length = 0
		if s[0:2] == '//':
			skipped = s.find('\n') + 1
			s = s[skipped:]

		start = s.find("<")
		end = s.find(">")
		try:
			name, length = s[start + 1:end].split(':')
			length = int(length)
		except ValueError:
			name = s[start + 1:end]

		value = s[end + 1:end + length + 1] if length else None
		return name, value, skipped + end + length + 1

	@staticmethod
	def parse_adif(stream):
		"""Consume ADIF tags out of @stream until an <EOR> tag is found.
		Returns a dictionary of the parsed key/value pairs, and the
		number of characters consumed from @stream."""

		tag = None
		off = 0
		d = {}
		while tag not in ["eor", "EOR"] and '<' in stream[off:]:
			tag, val, skip = QSO.next_adif_tag(stream[off:])
			d[tag.lower()] = val
			off += skip

		return d, off

	@staticmethod
	def get_qsos_from_file(path, prev=None, only=False):
		"""Returns a dictionary of QSOs keyed by callsign, containing
		all the QSO records that could be parsed from @path."""

		ret = prev or defaultdict(list)
		off = 0
		with open(path, 'r') as f:
			data = f.read().replace('\n', ' ')
			off = data.find("<eoh>") + 6
			valid = True

			while '<' in data[off:] and '>' in data[off:]:
				d, adj = QSO.parse_adif(data[off:])
				off += adj

				try:
					qso = QSO(d)
					if prev is not None and only \
							and qso.call not in ret:
						continue

					ret[qso.call].append(qso)
				except Exception as e:
					print("%s during parsing: '%s' (%s)" % (
					      e.__class__.__name__, str(e),
					      json.dumps(d)), file=sys.stderr)

		return ret

	@staticmethod
	def load_qsos_from_files(*args, ignore_invalid=True):
		"""Returns merged QSO records from multiple ADIF files, where a
		records with the same callsign and band are combined if their
		timestamps are within FUZZTIME of each other. If multiple files
		have different values for the same field, the last file in *args
		wins."""

		raw = None
		for filename in args:
			raw = QSO.get_qsos_from_file(filename, raw, True)

		qsos = []
		for call, qsolist in raw.items():
			for qcall in set(q.call for q in qsolist):
				fil = filter(lambda q: q.call == qcall, qsolist)
				o = list(sorted(fil, key=lambda q: (q.band,
							 q.dateon)))

				for i in range(len(o) - 1):
					a, b = o[i], o[i + 1]

					if b.dateon - a.dateon < QSO.FUZZTIME \
					   and a.band == b.band:
						o[i + 1] = a.merge(b)
						o[i] = None

				qsos.extend(filter(lambda x: x is not None, o))

		return qsos

	def merge(self, other):
		self.adifdict.update(other.adifdict)
		return self

	def parse_gridvector(self):
		src = self.my_gridsquare
		dst = self.gridsquare
		if not src or not dst:
			return

		# Truncate to 6-letter grid
		src, dst = src[:6].upper(), dst[:6].upper()

		# For 4x4, guess middle of gridsquare
		if len(src) == 4:
			src += "NN"
		if len(dst) == 4:
			dst += "NN"

		# A0A..R9X == 0..4319
		def v(s): return ord(s) - ord('A')
		srclon = v(src[0]) * (10 * 24) + int(src[2]) * 24 + v(src[4])
		srclat = v(src[1]) * (10 * 24) + int(src[3]) * 24 + v(src[5])
		dstlon = v(dst[0]) * (10 * 24) + int(dst[2]) * 24 + v(dst[4])
		dstlat = v(dst[1]) * (10 * 24) + int(dst[3]) * 24 + v(dst[5])

		# Normalize LON on [-pi, pi] and LAT on [-pi/2, pi/2]
		srclon = (srclon / 4320.0) * (2 * math.pi) - math.pi
		srclat = (srclat / 4320.0) * (math.pi) - (math.pi / 2)
		dstlon = (dstlon / 4320.0) * (2 * math.pi) - math.pi
		dstlat = (dstlat / 4320.0) * (math.pi) - (math.pi / 2)
		lmbda = dstlon - srclon

		# https://en.wikipedia.org/wiki/Great-circle_distance#Formulas
		rho = math.acos(math.sin(srclat) * math.sin(dstlat) +
				math.cos(srclat) * math.cos(dstlat) *
				math.cos(lmbda)
		)

		# https://en.wikipedia.org/wiki/Great-circle_navigation
		azrads = math.atan2(math.sin(lmbda),
				math.cos(srclat) * math.tan(dstlat) -
				math.sin(srclat) * math.cos(lmbda)
		)
		azdegs = (azrads / (math.pi * 2) * 360 + 360) % 360

		self.adifdict.update({
			"dist": int(abs(rho * 3959)) or 1, # Radius of Earth (miles)
			"azdegs": int(azdegs),
			"azrads": float(azrads),
			"cardinal": str(QSO.CARDINALS[math.floor(azdegs / 22.5)]),
		})

	def parse_date(self):
		d = {}

		# YYYYMMDD HHmmSS
		def adifdate_to_datetime(d, t):
			return datetime.datetime(tzinfo=datetime.timezone.utc,
				year=int(d[0:4]), month=int(d[4:6]),
				day=int(d[6:8]), hour=int(t[0:2]),
				minute=int(t[2:4]), second=int(t[4:6]))

		d["dateon"] = adifdate_to_datetime(
			self.adifdict["qso_date"],
			self.adifdict.get("time_on", "000000"),
		)
		d["epochon"] = d["dateon"].timestamp()

		if "qso_date_off" in self.adifdict:
			d["dateoff"] = adifdate_to_datetime(
				self.adifdict["qso_date_off"],
				self.adifdict.get("time_off", "000000"),
			)
			d["epochoff"] = d["dateoff"].timestamp()

		if "qslrdate" in self.adifdict:
			d["dateqsl"] = adifdate_to_datetime(
				self.adifdict["qslrdate"], "000000",
			)
			d["epochqsl"] = d["dateqsl"].timestamp()

		self.adifdict.update(d)

	def normalize_standard_fields(self):
		for fld, dfl in QSO.DEFAULT_FIELDS.items():
			self.adifdict[fld] = self.adifdict.get(fld, dfl)

		for fld, typ in QSO.PARSED_FIELDS.items():
			if fld not in self.adifdict:
				continue

			self.adifdict[fld] = typ(self.adifdict[fld])

	def __init__(self, adifdict):
		self.adifdict = copy.deepcopy(adifdict)
		self.normalize_standard_fields()
		self.parse_gridvector()
		self.parse_date()

	def __eq__(self, other):
		return abs(self.epochon - other.epochon) <= TIMEFUZZ \
		       and self.call == other.call and self.band == other.band

	def __bool__(self): return self.qsl_rcvd
	def __getattr__(self, name): return self.adifdict.get(name, None)
	def __getitem__(self, key): return self.adifdict.get(key, None)
	def __iter__(self): return iter(self.adifdict.items())
	def __contains__(self, item): return item in adifdict

class QSOGroup(object):
	def __init__(self, qsos, call=None, mode=None, grid=None):
		self.qsos = [q for q in qsos]
		if call:
			self.qsos = [q for q in self.qsos if q.station_callsign == call]
		if mode:
			self.qsos = [q for q in self.qsos if q.mode == mode]
		if grid:
			self.qsos = [q for q in self.qsos if q.my_gridsquare == grid]

		self.qsls = [qso for qso in self.qsos if qso.qsl_rcvd]

	def __iter__(self): return iter(self.qsos)
	def __len__(self): return len(self.qsos)

def do_list(args):
	qsos = QSOGroup(qsos=QSO.load_qsos_from_files(*args.files),
			call=args.call, mode=args.mode)
	keys = defaultdict(int)
	for qso in qsos:
		for key, _ in qso:
			keys[key] += 1

	# By default, just dump the list of possible fields
	if not args.fields:
		for key in sorted(keys.keys()):
			print("% 30s: %d/%d (%.3f%%)" % (key, keys[key],
			      len(qsos), keys[key] / len(qsos) * 100))
		return 1

	# Build all the columns/rows to print
	cols = [[str(q[x]) for x in args.fields] for q in
		sorted(qsos, key=lambda q: [q[x] for x in args.sorts])]

	# If a delimiter is specified, we don't do any padding
	if args.delim:
		print('\n'.join([args.delim.join(row) for row in cols]))
		return 0

	# Pretty table: pad all the columns so the longest value fits
	colwidths = [max(len(row[i]) for row in cols) for i in range(len(cols[0]))]
	for row in cols:
		fmt = " ".join([("%% %ds" % (x)) for x in colwidths])
		print(fmt % tuple(row))

	return 0

def export_adif(qsos):
	print(datetime.datetime.now(datetime.timezone.utc))
	print("<programid:8>adifutil")
	print("<eoh>", end="\n\n")

	for qso in qsos:
		for key, val in qso.adifdict.items():
			k = str(key)
			v = str(val)
			print("<%s:%d>%s" % (k, len(v), v))
		print("<eor>", end="\n\n")

	return 0

def export_json(qsos):
	types = [bool, int, float, str, type(None)]
	jsonout = [{str(k): (str(v) if type(v) not in types else v)
		    for k, v in iter(qso)} for qso in qsos]
	print(json.dumps(jsonout, indent="  "))
	return 0

def do_export(args):
	return {
		"adif": export_adif,
		"json": export_json,
	}[args.format](QSOGroup(qsos=QSO.load_qsos_from_files(*args.files)))

def parse_arguments():
	parser = argparse.ArgumentParser()
	parser.add_argument("files", type=str, nargs="+")
	subs = parser.add_subparsers(dest="action")

	listcmd = subs.add_parser("list")
	listcmd.add_argument("-d", dest="delim", type=str)
	listcmd.add_argument("-M", dest="mode", type=str, default=None)
	listcmd.add_argument("-C", dest="call", type=str, default=None)
	listcmd.add_argument("-F", dest="fields", type=str, nargs="*")
	listcmd.add_argument("-S", dest="sorts", type=str, nargs="*",
			     default=["dateon"])

	expcmd = subs.add_parser("export")
	expcmd = expcmd.add_argument("-F", dest="format", default="adif",
				     choices=["adif", "json"])

	return parser.parse_args()

def main():
	args = parse_arguments()
	return {
		"list": do_list,
		"export": do_export,
	}[args.action](args)

if __name__ == "__main__":
	sys.exit(main())
