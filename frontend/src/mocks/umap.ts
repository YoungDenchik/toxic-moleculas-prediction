import type { UmapPoint } from "../types";

export const mockUmap: UmapPoint[] = Array.from({ length: 200 }, (_, i) => ({
  id: `mock_${i}`,
  x: (Math.random() - 0.5) * 6,
  y: (Math.random() - 0.5) * 6,
  label: Math.random() > 0.6 ? 1 : 0,
  smiles: "CCN(CC)CCOc1ccc2nc(S(N)(=O)=O)sc2c1",
}));
