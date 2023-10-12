# from . import pyballcools
import os,sys
import os.path
import pandas as pd
from .pyballcools import (
	BAllCIndex,
	BAllC,
	IndexBallc,
	AllCToBallC,
	BAllCToAllC,
	ExtractCMeta,
	IndexCMeta
)
import pysam
import fire

import struct
from .bgzf import BgzfWriter,BgzfReader
dtype_func={
		'h': int, 'H': int,
		'i':int, 'I':int,'b': int, 'B':int,
		'L':int,'l':int,'q':int,'Q':int,
		'f':float,'d':float,
		's':str,'c':str
		}

class BAllCFile:
	def __init__(self, ballc_file, cmeta_file=None):
		"""
		Python warpper of BallCFile

		Parameters
		----------
		ballc_file: str
			path for ballc file (should be indexed)
		cmeta_file: str
			path for cmeta file ((should be indexed))
		"""
		self.bci = BAllCIndex(ballc_file)
		# self.tbi = tabix.open(cmeta_file)  if cmeta_file is not None else None
		self.tbi = pysam.TabixFile(cmeta_file) if cmeta_file is not None else None
		self.ballc=BAllC(ballc_file,"r")

	def __enter__(self):
		return self

	def __exit__(self, exc_type, exc_val, exc_tb):
		pass

	def header(self):
		if hasattr(self,'header_dict'):
			return self.header_dict
		self.header_dict={}
		attrs = ['version_minor', 'sc', 'assembly_text', 'l_assembly', 'header_text', 'l_text', 'refs', 'n_refs']
		for attr in attrs:
			# print(attr, self.ballc.header.__getattribute__(attr))
			self.header_dict[attr]=self.ballc.header.__getattribute__(attr)
		return self.header_dict

	def _fetch_with_cmeta(self, chrom, start, end):
		if chrom=="*":
			mciter = self.bci.QueryMcRecords_Iter("*")
		else:
			mciter = self.bci.QueryMcRecords_Iter(chrom, start, end)
		if mciter.HasNext():
			mciter.Next()
		while mciter.HasNext():
			rec = mciter.Next()
			try:
				record = next(self.tbi.fetch(rec.chrom,rec.pos-1,rec.pos,)).split()
				# record = next(self.tbi.query(rec.chrom,rec.pos-1,rec.pos))
				*_, strand, context = record
				yield(rec.chrom,rec.pos,strand, context, rec.mc,rec.cov, )
			except:
				pass
				# print(f"No meta data found for {rec.chrom}:{rec.pos-1}-{rec.pos}")

	def _fetch(self, chrom, start, end):
		if chrom=="*":
			mciter = self.bci.QueryMcRecords_Iter("*")
		else:
			mciter = self.bci.QueryMcRecords_Iter(chrom, start, end)
		if mciter.HasNext():
			mciter.Next()
		while mciter.HasNext():
			rec = mciter.Next()
			yield(rec.chrom, rec.pos, rec.mc, rec.cov, )

	def fetch(self, chrom, start, end):
		"""
		Fetch region from ballc.
		Parameters
		----------
		chrom: str
			chromsome name,  if chrom=="*", fetch all records.
		start: int
			start position, could be None if chrom=='*'
		end: int
			end position, could be None if chrom=='*'

		Returns
		-------
		a generator containing the records.
		"""
		if self.tbi is None:
			return self._fetch(chrom, start, end)
		else:
			return self._fetch_with_cmeta(chrom, start, end)

	def fetch_line(self, chrom, start, end):
		if self.tbi is None:
			mciter = self._fetch(chrom, start, end)
			for rec in mciter:
				yield('{}\t{}\t{}\t{}'.format(*rec))
		else:
			mciter = self._fetch_with_cmeta(chrom, start, end)
			for rec in mciter:
				yield('{}\t{}\t{}\t{}\t{}\t{}'.format(*rec))

	def to_allc(self,allc_path):
		allc_path=os.path.abspath(os.path.expanduser(allc_path))
		f=open(allc_path,'w')
		for line in self.fetch_line("*",None,None):
			f.write(line+"\n")
		f.close()
		return allc_path

def Ballc2Allc(ballc_path,cmeta_path,
			   allc_path,warn_mismatch=True,
			   err_mismatch=True,skip_mismatch=True,
			   c_context="*"):
	"""
	Convert ballc file into allc path.

	Parameters
	----------
	ballc_path: str
		input ballc path, should be indexed
	cmeta_path:str
		path
	allc_path: str
		output allc file
	warn_mismatch: bool
		warn_mismatch
	err_mismatch: bool
		err_mismatch
	skip_mismatch: bool
		skip_mismatch
	c_context: str
		c_context

	Returns
	-------

	"""
	# bf = BAllCFile(ballc_path, cmeta_path)
	# allc_path=bf.to_allc(allc_path)
	BAllCToAllC(ballc_path, cmeta_path, allc_path,
				 warn_mismatch, err_mismatch, skip_mismatch,
				 c_context)
	return allc_path

def index_ballc(ballc_path):
	IndexBallc(ballc_path)

def Allc2Ballc(allc_path,ballc_path,chrom_size_path=None,
				assembly_text="",header_text="",sc=True):
	"""
	Convert allc file into ballc file.

	Parameters
	----------
	allc_path: str
		input allc file path.
	ballc_path: str
		output ballc path, will be indexed automatically.
	chrom_size_path: str
		path
	assembly_text: str
		text to be added
	header_text: str
		text to be added
	sc: bool
		whether single cell file?

	Returns
	-------

	"""
	AllCToBallC(allc_path,ballc_path,chrom_size_path,
				assembly_text,header_text,sc)
	index_ballc(ballc_path)
	return ballc_path

def extractC(fasta_path,cmeta_path):
	"""
	Extract all C position from fasta file.

	Parameters
	----------
	fasta_path: str
		path for fasta file
	cmeta_path: str
		path for the output cmeta file.

	Returns
	-------

	"""
	ExtractCMeta(fasta_path, cmeta_path)
	IndexCMeta(cmeta_path)
	return cmeta_path

def header(ballc_path,cmeta_path):
	"""
	Print ballc file header.

	Parameters
	----------
	ballc_path: str
		path for input ballc path
	cmeta_path: str
		path for input cmeta_path

	Returns
	-------

	"""
	bf = BAllCFile(ballc_path, cmeta_path)
	header_dict=bf.header()
	for key in header_dict:
		value=header_dict[key]
		print(f"{key}: {value}")

def query(ballc_path,cmeta_path=None,
		  chrom="*", start=None, end=None):
	"""
	Query ballc file with or without cmeta index.

	Parameters
	----------
	ballc_path: str
		path for ballc file.
	cmeta_path: str
		path for cmeta file
	chrom: str
		chromosome, "*" to query all records.
	start: int
		start position, if chrom=="*", start can be ignored.
	end: int
		start position, if chrom=="*", start can be ignored.

	Returns
	-------

	"""
	bf = BAllCFile(ballc_path, cmeta_path)
	for line in bf.fetch_line(chrom, start, end):
		try:
			sys.stdout.write(line+'\n')
		except:
			sys.stdout.close()
			break

def main():
	fire.core.Display = lambda lines, out: print(*lines, file=out)
	fire.Fire({
		"cmeta":extractC,
		"b2a":Ballc2Allc,
		"a2b":Allc2Ballc,
		"header":header,
		"query": query
	})

def open_ballc(ballc_path,lib='bgzf'):
	if lib=='bgzf':
		from Bio import bgzf
		return bgzf.open(ballc_path, 'rb')
	elif lib=='bgzip':
		import bgzip
		return bgzip.BGZipReader(open(ballc_path,'rb'))
	elif lib=='gzip':
		import gzip
		return gzip.open(ballc_path,'rb')
	else:
		return open(ballc_path,'rb')
	
def fread(f,size,fmt):
	return struct.unpack(fmt, f.read(size))

def read_header(f):
	D={}
	if isinstance(f,str):
		f=open_ballc(f)
	D['magic']=fread(f,6,'6s')[0].decode()
	D['version'],D['version_minor'],D['sc']=fread(f,3,'3B')
	D['l_assembly']=fread(f,4,'I')[0]
	if D['l_assembly'] > 0:
		text=fread(f,D['l_assembly'],f"{D['l_assembly']}s")[0].decode()
	else:
		text=''
	D['assembly_text']=text
	D['l_text'] = fread(f, 4, 'I')[0]
	text=fread(f, D['l_text'], f"{D['l_text']}s")[0].decode()
	D['text'] = text
	D['n_refs']=fread(f, 2, 'H')[0]
	D['refs'] = []
	D['ref_length'] = {}
	for i in range(D['n_refs']):
		l_name=fread(f, 4, 'I')[0]
		ref_name=fread(f, l_name, f"{l_name}s")[0].decode()
		ref_len=fread(f, 4, 'I')[0]
		D['refs'].append(ref_name)
		D['ref_length'][ref_name]=ref_len
	D['offset']=f.tell()
	return D

def get_ballc_records(f,offset,refs,sc=1,chunk_size=50000):
	f.seek(offset)
	if sc:
		fmts='<I3H'
	else:
		fmts='<IH2B'
	size=struct.calcsize(fmts) * chunk_size
	# pos=fread(f, 4, 'I')[0]
	# ref_id=fread(f, 2, 'H')[0]
	# mc=fread(f, 2, 'H')[0] if sc else fread(f, 1, 'B')[0]
	# cov=fread(f, 2, 'H')[0] if sc else fread(f, 1, 'B')[0]
	flag=1
	while flag:
		chunk=struct.iter_unpack(fmts,f.read(size))
		while True:
			try:
				r=chunk.__next__()
				yield (refs[r[1]], r[0], r[2], r[3])
			except StopIteration:
				flag=0
				break
	
def read_ballc(ballc_path):
	f=open_ballc(ballc_path)
	D=read_header(f) #f._buffer, f.block_raw_length, f._block_start_offset
	# print(D)
	records=get_ballc_records(f,D['offset'],D['refs'],sc=1)
	for record in records:
		print(record)
	
def read_bci(bci_path,chunk_size=5000):
	D={}
	f=open_ballc(bci_path) #magic (1 byte); size (8,Q), (ref_id (2bytes int, H), bin (4 bytes int, I), chunk_start,chunk_end (uint64_t, 8bytes,Q))
	D['magic']=fread(f,8,'<8s')[0].decode()
	D['size']=fread(f,8,'<Q')[0]
	D['chunks']=[]
	fmts="<HIQQ"
	fmt_size = struct.calcsize(fmts) * chunk_size
	flag = 1
	while flag:
		chunk = struct.iter_unpack(fmts, f.read(fmt_size))
		while True:
			try:
				r = chunk.__next__()
				D['chunks'].append(r)
			except StopIteration:
				flag = 0
				break
	f.close() #[(0, 4864, 1432, 13232), (0, 4865, 13242, 28292), (0, 4866, 28302, 44132)]
	return D
	
def stdin_generator(fmts,sep='\t'):
	line=sys.stdin.readline()
	values=line.split(sep)
	if type(fmts)==str and len(fmts)==1:
		fmts=fmts*len(values)
	functions=[dtype_func[f] for f in fmts]
	yield [func(v) for v,func in zip(values,functions)]
	while line:
		line = sys.stdin.readline()
		values = line.split(sep)
		yield struct.pack(fmts,*[func(v) for v, func in zip(values, functions)])
		
def df_generator(df,fmts):
	if isinstance(df,pd.DataFrame):
		values=df.values
	else:
		values=df
	if type(fmts)==str and len(fmts)==1:
		fmts=fmts*len(values[0])
	for value in values:
		yield struct.pack(fmts,*list(value))
	
def pack(Input, fmts='H',sep='\t',Output=None,usecols=None,
		 skiprows=0,chunksize=5000):
	"""
	Pack dataframe or file into binary bgzf. pack("10_merged_allc.tsv.gz",
		fmts='HH',Output="test.bgzf",usecols=[4,5])
	Parameters
	----------
	Input :
	fmts :
	sep :
	Output :
	usecols :
	skiprows :

	Returns
	-------

	"""
	if not Output is None:
		Output=os.path.expanduser(Output)
	else:
		Output="out.bgzf"
	writer=BgzfWriter(Output,mode="w",format=fmts)
	if Input is None: #input from stdin
		# data = stdin_generator(fmts,sep)
		stdin=sys.stdin.buffer
		while True:
			data = stdin.read(1 << 16) #65536=2**16
			writer.write(data)
			if not data:
				break
	else:
		if isinstance(Input,str):
			input_path=os.path.expanduser(Input)
			dfs=pd.read_csv(input_path,sep=sep,usecols=usecols,
						   header=None,skiprows=skiprows,chunksize=chunksize)
			if type(fmts)==str and len(fmts)==1:
				fmts=fmts*dfs[0].shape[1]
			for df in dfs:
				writer.write(df.apply(lambda x: struct.pack(fmts, *x.tolist()), axis=1).sum())
		else: #Input is a dataframe
			df=Input
			if type(fmts)==str and len(fmts)==1:
				fmts=fmts*df.shape[1]
			writer.write(df.apply(lambda x: struct.pack(fmts, *x.tolist()), axis=1).sum())
	writer.close()
	
def view(Input="test.bgzf"):
	reader=BgzfReader(Input,'rb')
	while reader._block_raw_length > 0:
		for r in struct.iter_unpack(reader.format,reader._buffer):
			yield r
		reader._load_block()
	reader.close()

def create_index(Input):
	index=[]
	for block in bgzf.BgzfBlocks(open(Input, 'rb')):
		index.append(block)
	return index

if __name__=="__main__""" :
	main()