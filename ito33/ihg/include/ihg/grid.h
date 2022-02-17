#ifndef GRID_H
#define GRID_H




namespace ito33
{

namespace ihg
{

enum SpecialPointType
{
  eRefineLeft,
  eRefineRight,
  eRefineBoth,
  eRefineNone
};


class Grid
{
protected:
  double* mpdGrid;
  double mdGrowth;
  int miNumPoints;

  double* m_pdSpecialPoints;
  SpecialPointType* m_pSpecialPointTypes;
  size_t m_nNbSpecialPoints;
  double* workArray;
  void Construct(size_t nNbPts);
  void Equidistribute(size_t nNbPts, size_t nMaxIteration, double* pdGrid, double dLeft, double dRight);
  void GetDensities(double* pdX, double* pdWeights, size_t nNbX);
  bool TridiagSolve(double* a, double* b, double* c, double* r, int n, double* x);



  void ComputeWidths(int iNumPoints, double* dPoints, SpecialPointType* eTypes, double* dWidthLeft, double* dWidthRight);

public:

  Grid(int iNpts, double* pdSpecialPoints, SpecialPointType* pSpecialPointTypes, int iNSpecialPoints, size_t nRefine);
  ~Grid();

  int GetNumPoints(void);
  double* GetGrid(void);
  void UniformRefine(void);
  
};

} // namespace ihg

} // namespace ito33

#endif

