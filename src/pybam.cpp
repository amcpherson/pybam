#include <iostream>
#include <iomanip>
#include <vector>
#include <queue>

#include <boost/foreach.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_utility.hpp>
#include <boost/graph/incremental_components.hpp>
#include <boost/graph/breadth_first_search.hpp>
#include <boost/pending/disjoint_sets.hpp>
#include <boost/graph/strong_components.hpp>
#include <boost/lexical_cast.hpp>

#include <boost/python/class.hpp>
#include <boost/python/module.hpp>
#include <boost/python/def.hpp>
#include <boost/python/list.hpp>
#include <boost/python/tuple.hpp>

#include "bamtools/src/api/BamReader.h"
#include "bamtools/src/utils/bamtools_pileup_engine.h"
#include "bamtools/src/utils/bamtools_fasta.h"
#include "bamtools/src/utils/bamtools_utilities.h"

using namespace std;
using namespace boost;
using namespace BamTools;


python::tuple CreatePileupTuple(const PileupPosition& pileupData)
{
	int ntData[5][6] = {{0}};
	int ambiguous = 0;
	int insertionCount = 0;
	int deletionCount = 0;
	
	for (vector<PileupAlignment>::const_iterator pileupIter = pileupData.PileupAlignments.begin(); pileupIter != pileupData.PileupAlignments.end(); ++pileupIter)
	{
		const PileupAlignment& pa = (*pileupIter);
		const BamAlignment& ba = pa.Alignment;
		
		// adjacent insertions and deletions
		for (std::vector<CigarOp>::const_iterator opIter = ba.CigarData.begin(); opIter != ba.CigarData.end(); opIter++)
		{
			if (opIter->Type == 'I')
			{
				insertionCount++;
			}
			else if (opIter->Type == 'D')
			{
				deletionCount++;
			}
		}
		
		if (pa.IsCurrentDeletion)
		{
			continue;
		}
		
		char base = toupper(ba.QueryBases.at(pa.PositionInAlignment));
		
		if (base == 'N')
		{
			ambiguous++;
			continue;
		}
		
		int baseIdx;
		switch (base)
		{
			case 'A': baseIdx = 0; break;
			case 'C': baseIdx = 1; break;
			case 'G': baseIdx = 2; break;
			case 'T': baseIdx = 3; break;
			default: throw runtime_error("unrecognized base " + string(1, base));
		}
		
		// count
		ntData[baseIdx][0]++;
		ntData[4][0]++;
		
		// quality
		ntData[baseIdx][1] += ba.Qualities.at(pa.PositionInAlignment);
		ntData[4][1] += ba.Qualities.at(pa.PositionInAlignment);
		
		// mapping quality
		ntData[baseIdx][2] += ba.MapQuality;
		ntData[4][2] += ba.MapQuality;
		
		// distance
		ntData[baseIdx][3] += (ba.IsReverseStrand()) ? ba.Length - pa.PositionInAlignment - 1 : pa.PositionInAlignment;
		ntData[4][3] += (ba.IsReverseStrand()) ? ba.Length - pa.PositionInAlignment - 1 : pa.PositionInAlignment;
		
		// direction
		ntData[baseIdx][3] += (ba.IsReverseStrand()) ? 1 : 0;
		ntData[4][3] += (ba.IsReverseStrand()) ? 1 : 0;
	}
	
	// Identify major base
	int majorBaseIdx = 0;
	for (int baseIdx = 0; baseIdx < 4; baseIdx++)
	{
		if (ntData[baseIdx][0] > ntData[majorBaseIdx][0])
		{
			majorBaseIdx = baseIdx;
		}
	}
	
	// Identify minor base, initialize to base that is not the major base
	int minorBaseIdx = (majorBaseIdx + 1) % 4;
	for (int baseIdx = 0; baseIdx < 4; baseIdx++)
	{
		if (ntData[baseIdx][0] > ntData[minorBaseIdx][0] && baseIdx != majorBaseIdx)
		{
			minorBaseIdx = baseIdx;
		}
	}
	
	// Calculate entropy
	double depth = (double)pileupData.PileupAlignments.size(); // NOTE: COPIED FROM JEFF, THIS MAY BE A WRONG
	double entropy = 0.0;
	for (int baseIdx = 0; baseIdx < 4; baseIdx++)
	{
		double probability = (double)ntData[baseIdx][0] / depth;
		if (probability != 0)
		{
			entropy -= (log(probability) * probability);
		}
	}
	
	// Interface is 1-based, bamtools is 0-based
	int position = pileupData.Position + 1;
	
	return python::make_tuple(position,
							  python::make_tuple(ntData[0][0], ntData[0][1], ntData[0][2], ntData[0][3], ntData[0][4], ntData[0][5]),
							  python::make_tuple(ntData[1][0], ntData[1][1], ntData[1][2], ntData[1][3], ntData[1][4], ntData[1][5]),
							  python::make_tuple(ntData[2][0], ntData[2][1], ntData[2][2], ntData[2][3], ntData[2][4], ntData[2][5]),
							  python::make_tuple(ntData[3][0], ntData[3][1], ntData[3][2], ntData[3][3], ntData[3][4], ntData[3][5]),
							  python::make_tuple(ntData[4][0], ntData[4][1], ntData[4][2], ntData[4][3], ntData[4][4], ntData[4][5]),
							  majorBaseIdx,
							  minorBaseIdx,
							  ambiguous,
							  insertionCount,
							  entropy,
							  deletionCount,
							  pileupData.RefId
							  );
}


struct PileupQueue : PileupVisitor
{
	PileupQueue() : StartRefId(-1), StartPosition(-1) {}
	
	void Visit(const PileupPosition& pileupData)
	{
		// Reset if we ended up on the wrong chromosome
		if (StartRefId >= 0 && pileupData.RefId != StartRefId)
		{
			StartRefId = -1;
			StartPosition = -1;
		}
		
		// Dont store tuples before the start position
		if (pileupData.Position < StartPosition)
		{
			return;
		}
		
		Pileups.push(CreatePileupTuple(pileupData));
		
		// Reset start refid/position to have no further effect
		StartRefId = -1;
		StartPosition = -1;
	}
	
	void Clear()
	{
		std::queue<python::tuple> empty;
		swap(Pileups, empty);
	}
	
	std::queue<python::tuple> Pileups;
	int StartRefId;
	int StartPosition;
};

class PyPileup
{
public:
	PyPileup() : m_PileupEngine(0), m_PileupQueue(0)
	{
	}
	
	~PyPileup()
	{
		delete m_PileupEngine;
		delete m_PileupQueue;
	}
	
	void Open(const string& bamFilename)
	{
		if (!m_BamReader.Open(bamFilename))
		{
			throw runtime_error("unable to open bam file " + bamFilename);
		}
		
		if (!m_BamReader.LocateIndex())
		{
			throw runtime_error("unable to open index for bam file " + bamFilename);
		}
		
		RefNames = python::list();
		for (RefVector::const_iterator refDataIter = m_BamReader.GetReferenceData().begin(); refDataIter != m_BamReader.GetReferenceData().end(); refDataIter++)
		{
			RefNames.append(refDataIter->RefName);
		}
		
		RestartPileupEngine();
	}
	
	void Rewind()
	{
		m_BamReader.Rewind();
		
		RestartPileupEngine();
	}
	
	void JumpRef(const string& refName)
	{
		int refId = m_BamReader.GetReferenceID(refName);
		
		if (refId < 0)
		{
			throw runtime_error("invalid ref name " + refName);
		}
		
		m_BamReader.Jump(refId);
		
		RestartPileupEngine();
	}
	
	void JumpRefPosition(const string& refName, int position)
	{
		// Interface is 1-based, bamtools is 0-based
		position -= 1;
		
		int refId = m_BamReader.GetReferenceID(refName);
		
		if (refId < 0)
		{
			throw runtime_error("invalid ref name " + refName);
		}
		
		m_BamReader.Jump(refId, position);
		
		RestartPileupEngine();
		
		m_PileupQueue->StartRefId = refId;
		m_PileupQueue->StartPosition = position;
	}
	
	python::object Next()
	{
		if (m_PileupEngine == 0)
		{
			throw runtime_error("next called before open");
		}
		
		BamAlignment al;
		while (m_BamReader.GetNextAlignment(al))
		{
			m_PileupEngine->AddAlignment(al);
			
			if (!m_PileupQueue->Pileups.empty())
			{
				return PopPileup();
			}
		}
		m_PileupEngine->Flush();
		
		if (!m_PileupQueue->Pileups.empty())
		{
			return PopPileup();
		}
		
		return python::object();
	}
	
	python::list RefNames;
	
private:
	void RestartPileupEngine()
	{
		delete m_PileupEngine;
		m_PileupEngine = new PileupEngine();
		
		delete m_PileupQueue;
		m_PileupQueue = new PileupQueue();
		
		m_PileupEngine->AddVisitor(m_PileupQueue);
	}
	
	python::object PopPileup()
	{
		python::tuple tpl;
		swap(tpl, m_PileupQueue->Pileups.front());
		m_PileupQueue->Pileups.pop();
		
		return tpl;
	}
	
	BamReader m_BamReader;
	
	PileupEngine* m_PileupEngine;
	PileupQueue* m_PileupQueue;
};

class PyFasta
{
public:
	PyFasta() : m_IsOpen(false)
	{
	}
	
	~PyFasta()
	{
		if (m_IsOpen)
		{
			m_Fasta.Close();
		}
	}
	
	void Open(const string& fastaFilename)
	{
		if (!Utilities::FileExists(fastaFilename))
		{
			throw runtime_error("invalid fasta file " + fastaFilename);
		}
		
		string indexFilename = fastaFilename + ".fai";
		
		if (!Utilities::FileExists(indexFilename))
		{
			throw runtime_error("index file " + indexFilename + " not found");
		}
		
		if (!m_Fasta.Open(fastaFilename, indexFilename))
		{
			throw runtime_error("unable to open fasta file " + fastaFilename);
		}
		
		vector<string> referenceNames = m_Fasta.GetReferenceNames();
		for (int refId = 0; refId < referenceNames.size(); refId++)
		{
			m_RefNameId[referenceNames[refId]] = refId;
		}
		
		m_RefLengths = m_Fasta.GetReferenceLengths();
		
		m_IsOpen = true;
	}
	
	python::object GetPosition(const string& refName, int position)
	{
		// Interface is 1-based, bamtools is 0-based
		position -= 1;
		
		if (!m_IsOpen)
		{
			throw runtime_error("get called before open");
		}
		
		unordered_map<string,int>::const_iterator refNameIdIter = m_RefNameId.find(refName);
		if (refNameIdIter == m_RefNameId.end())
		{
			throw runtime_error("unknown ref name " + refName);
		}
		int refId = refNameIdIter->second;
		
		char referenceBase = 'N';
		if (!m_Fasta.GetBase(refId, position, referenceBase))
		{
			throw runtime_error("unable to get base at " + refName + ":" + lexical_cast<string>(position));
		}
		
		return python::make_tuple(referenceBase, 0, 0, 0, 0);
	}
	
private:
	Fasta m_Fasta;
	bool m_IsOpen;
	unordered_map<string,int> m_RefNameId;
	vector<int> m_RefLengths;
};

BOOST_PYTHON_MODULE(newpybam)
{
	using namespace python;
	
	class_<PyPileup>("pileup", init<>())
		.def_readonly("refnames", &PyPileup::RefNames)
		.def("open", &PyPileup::Open)
		.def("rewind", &PyPileup::Rewind)
		.def("jump", &PyPileup::JumpRef)
		.def("jump", &PyPileup::JumpRefPosition)
		.def("next", &PyPileup::Next)
	;
	
	class_<PyFasta>("fasta", init<>())
		.def("open", &PyFasta::Open)
		.def("get", &PyFasta::GetPosition)
	;
}

