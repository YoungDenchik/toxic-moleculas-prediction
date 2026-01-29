import Plot from "react-plotly.js";
import type { UmapPoint } from "../types";

export function UmapPlot(props: {
  points: UmapPoint[];
  onPick?: (p: UmapPoint) => void;
}) {
  const xs = props.points.map((p) => p.x);
  const ys = props.points.map((p) => p.y);

  // кольори по label (0/1)
  const colors = props.points.map((p) => (p.label === 1 ? 1 : 0));
  const texts = props.points.map((p) => `${p.id}${p.smiles ? `<br>${p.smiles}` : ""}`);

  return (
    <div style={{ border: "1px solid #e2e8f0", borderRadius: 16, padding: 10 }}>
      <div style={{ fontWeight: 800, marginBottom: 8 }}>UMAP</div>

      <Plot
        data={[
          {
            type: "scattergl",        // важливо: тягне багато точок краще
            mode: "markers",
            x: xs,
            y: ys,
            text: texts,
            hoverinfo: "text",
            marker: { size: 7, color: colors },
          } as any,
        ]}
        layout={{
          height: 420,
          margin: { l: 40, r: 20, t: 10, b: 40 },
          xaxis: { title: "UMAP-1" },
          yaxis: { title: "UMAP-2" },
        }}
        config={{ displayModeBar: true, responsive: true }}
        onClick={(e: any) => {
          const idx = (e?.points?.[0] as any)?.pointIndex;
          if (idx != null && props.onPick) props.onPick(props.points[idx]);
        }}
      />
      <div style={{ color: "#64748b", fontSize: 12 }}>
        Click a point to select it.
      </div>
    </div>
  );
}
