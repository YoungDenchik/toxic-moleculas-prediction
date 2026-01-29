import Plot from "react-plotly.js";
import type { PlotMouseEvent } from "plotly.js";
import type { UmapPoint, ChemicalSpace } from "../types";

interface UmapPlotProps {
  points: UmapPoint[];
  onPick?: (p: UmapPoint) => void;
  userMolecule?: {
    position: ChemicalSpace;
    smiles: string;
    isToxic: boolean;
  };
}

type PlotData = React.ComponentProps<typeof Plot>["data"][number];

export function UmapPlot({ points, onPick, userMolecule }: UmapPlotProps) {
  // Separate toxic and non-toxic points
  const toxicPoints = points.filter((p) => p.label === 1);
  const nonToxicPoints = points.filter((p) => p.label !== 1);

  const traces: PlotData[] = [
    // Non-toxic points (label 0)
    {
      type: "scattergl",
      mode: "markers",
      name: "Non-toxic",
      x: nonToxicPoints.map((p) => p.x),
      y: nonToxicPoints.map((p) => p.y),
      text: nonToxicPoints.map((p) => `${p.id}${p.smiles ? `<br>${p.smiles}` : ""}`),
      hoverinfo: "text",
      marker: {
        size: 6,
        color: "#22c55e",
        opacity: 0.6,
      },
    } as PlotData,
    // Toxic points (label 1)
    {
      type: "scattergl",
      mode: "markers",
      name: "Toxic",
      x: toxicPoints.map((p) => p.x),
      y: toxicPoints.map((p) => p.y),
      text: toxicPoints.map((p) => `${p.id}${p.smiles ? `<br>${p.smiles}` : ""}`),
      hoverinfo: "text",
      marker: {
        size: 6,
        color: "#ef4444",
        opacity: 0.6,
      },
    } as PlotData,
  ];

  // Add user's molecule as a highlighted point
  if (userMolecule) {
    traces.push({
      type: "scatter",
      mode: "markers+text",
      name: "Your molecule",
      x: [userMolecule.position.x],
      y: [userMolecule.position.y],
      text: ["Your molecule"],
      textposition: "top center",
      textfont: {
        size: 12,
        color: "#1d4ed8",
      },
      hoverinfo: "text",
      hovertext: [`Your molecule<br>${userMolecule.smiles}`],
      marker: {
        size: 16,
        color: userMolecule.isToxic ? "#dc2626" : "#16a34a",
        symbol: "star",
        line: {
          color: "#1d4ed8",
          width: 3,
        },
      },
    } as PlotData);
  }

  const handleClick = (e: PlotMouseEvent) => {
    const pointData = e.points?.[0];
    if (pointData && pointData.curveNumber < 2) {
      // Only for background points
      const idx = pointData.pointIndex;
      const targetPoints = pointData.curveNumber === 0 ? nonToxicPoints : toxicPoints;
      if (idx != null && onPick && targetPoints[idx]) {
        onPick(targetPoints[idx]);
      }
    }
  };

  return (
    <div style={{ border: "1px solid #e2e8f0", borderRadius: 16, padding: 12, background: "#fff" }}>
      <div style={{ fontWeight: 800, fontSize: 14, marginBottom: 8 }}>Chemical Space (UMAP)</div>

      <Plot
        data={traces}
        layout={{
          height: 400,
          margin: { l: 50, r: 20, t: 10, b: 50 },
          xaxis: {
            title: "UMAP-1",
            gridcolor: "#f1f5f9",
            zerolinecolor: "#e2e8f0",
          },
          yaxis: {
            title: "UMAP-2",
            gridcolor: "#f1f5f9",
            zerolinecolor: "#e2e8f0",
          },
          legend: {
            x: 0,
            y: 1.15,
            orientation: "h",
            bgcolor: "rgba(255,255,255,0.8)",
          },
          plot_bgcolor: "#fafafa",
          paper_bgcolor: "#fff",
          hovermode: "closest",
        }}
        config={{
          displayModeBar: true,
          responsive: true,
          displaylogo: false,
          modeBarButtonsToRemove: ["lasso2d", "select2d"],
        }}
        onClick={handleClick}
        style={{ width: "100%" }}
      />

      <div style={{ display: "flex", gap: 16, marginTop: 8, fontSize: 12, color: "#64748b" }}>
        <div style={{ display: "flex", alignItems: "center", gap: 4 }}>
          <span style={{ width: 10, height: 10, borderRadius: "50%", background: "#22c55e" }}></span>
          Non-toxic compounds
        </div>
        <div style={{ display: "flex", alignItems: "center", gap: 4 }}>
          <span style={{ width: 10, height: 10, borderRadius: "50%", background: "#ef4444" }}></span>
          Toxic compounds
        </div>
        {userMolecule && (
          <div style={{ display: "flex", alignItems: "center", gap: 4 }}>
            <span style={{ fontSize: 14 }}>⭐</span>
            Your molecule
          </div>
        )}
      </div>
    </div>
  );
}
