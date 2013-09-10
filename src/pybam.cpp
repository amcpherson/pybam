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

#include <boost/python/class.hpp>
#include <boost/python/module.hpp>
#include <boost/python/def.hpp>
#include <boost/python/list.hpp>
#include <boost/python/manage_new_object.hpp>

#include "bamtools/src/api/BamReader.h"
#include "bamtools/src/utils/bamtools_pileup_engine.h"

using namespace std;
using namespace boost;
using namespace BamTools;


struct PileupQueue : PileupVisitor
{
	void Visit(const PileupPosition& pileupData)
	{
		Pileups.push(pileupData);
	}
	
	void Clear()
	{
		std::queue<PileupPosition> empty;
		swap(Pileups, empty);
	}
	
	std::queue<PileupPosition> Pileups;
};

struct PileupInfo
{
	int position;
	string bases;
};

class Pileup
{
public:
	Pileup()
	{
		m_PileupEngine.AddVisitor(&m_PileupQueue);
	}
	
	void Open(const string& bamFilename)
	{
		if (!m_BamReader.Open(bamFilename))
		{
			throw runtime_error("Unable to open bam file " + bamFilename);
		}
		
		m_PileupEngine.Flush();
		m_PileupQueue.Clear();
	}
	
	PileupInfo* Next()
	{
		BamAlignment al;
		while (m_BamReader.GetNextAlignment(al) && m_PileupQueue.Pileups.empty())
		{
			m_PileupEngine.AddAlignment(al);
		}
		
		if (!m_PileupQueue.Pileups.empty())
		{
			const PileupPosition& pileupData = m_PileupQueue.Pileups.front();
			
			PileupInfo* pileupInfo = new PileupInfo();
			pileupInfo->position = pileupData.Position;
			
			vector<PileupAlignment>::const_iterator pileupIter = pileupData.PileupAlignments.begin();
			vector<PileupAlignment>::const_iterator pileupEnd  = pileupData.PileupAlignments.end();
			for (vector<PileupAlignment>::const_iterator pileupIter = pileupData.PileupAlignments.begin(); pileupIter != pileupData.PileupAlignments.end(); ++pileupIter)
			{
				const PileupAlignment& pa = (*pileupIter);
				const BamAlignment& ba = pa.Alignment;
				
				if ( !pa.IsCurrentDeletion )
				{
					char base = ba.QueryBases.at(pa.PositionInAlignment);
					pileupInfo->bases += base;
				}
			}
			
			m_PileupQueue.Pileups.pop();
			
			return pileupInfo;
		}
		
		return 0;
	}
	
private:
	BamReader m_BamReader;
	PileupEngine m_PileupEngine;
	PileupQueue m_PileupQueue;
};


BOOST_PYTHON_MODULE(pybam)
{
	using namespace python;
	
	class_<Pileup>("pileup", init<>())
		.def("open", &Pileup::Open)
		.def("next", &Pileup::Next, return_value_policy<manage_new_object>())
	;
	
	class_<PileupInfo>("pileupinfo", init<>())
		.def_readonly("position", &PileupInfo::position)
		.def_readonly("bases", &PileupInfo::bases)
	;
}

