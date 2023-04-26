#ifndef _DEL_TRANGULER_H_
#define _DEL_TRANGULER_H_


#include "subdivision2d.h"

class DelTriangulater
{
public:
	DelTriangulater();
	~DelTriangulater();

	MRESULT doinit(MLong lWidth, MLong lHeight);
	MRESULT uninit();
	MRESULT setTriangulePoints(std::vector<MPOINTF> &points);
	MRESULT triVertexToIndex(std::vector<MPOINTF> &points, std::vector<MUInt32> &vertex);
    MInt32  getTriangleNum() {return m_triNum;};
private:
	Subdiv2D*		m_subdiv;
    MInt32			m_triNum;
};
#endif // !_DEL_TRANGULE_H_

