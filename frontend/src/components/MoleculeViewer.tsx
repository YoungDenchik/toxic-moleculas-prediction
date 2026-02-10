import { useEffect, useState } from "react";

interface MoleculeViewerProps {
  smiles: string;
  width?: number;
  height?: number;
}

export function MoleculeViewer({ smiles, width = 300, height = 200 }: MoleculeViewerProps) {
  const [svgContent, setSvgContent] = useState<string>("");
  const [status, setStatus] = useState<"idle" | "loading" | "ready" | "error">("idle");
  const [errorMsg, setErrorMsg] = useState<string>("");

  useEffect(() => {
    if (!smiles) {
      setSvgContent("");
      setStatus("idle");
      return;
    }

    const cleanSmiles = smiles.trim();
    if (!cleanSmiles) {
      setStatus("error");
      setErrorMsg("Empty SMILES");
      return;
    }

    setStatus("loading");
    setErrorMsg("");

    // Use PubChem's REST API to get SVG image
    const encodedSmiles = encodeURIComponent(cleanSmiles);
    const imgUrl = `https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/${encodedSmiles}/PNG?image_size=${width}x${height}`;

    // Test if the image loads
    const img = new Image();
    img.onload = () => {
      setStatus("ready");
      setSvgContent(imgUrl);
    };
    img.onerror = () => {
      console.error("[MoleculeViewer] PubChem failed for:", cleanSmiles);
      setStatus("error");
      setErrorMsg("Could not render molecule");
    };
    img.src = imgUrl;
  }, [smiles, width, height]);

  if (!smiles) {
    return (
      <div
        style={{
          width,
          height,
          display: "flex",
          alignItems: "center",
          justifyContent: "center",
          border: "2px dashed #cbd5e1",
          borderRadius: 12,
          color: "#64748b",
          fontSize: 14,
        }}
      >
        Enter a SMILES to visualize
      </div>
    );
  }

  if (status === "error") {
    return (
      <div
        style={{
          width,
          height,
          display: "flex",
          alignItems: "center",
          justifyContent: "center",
          border: "1px solid #fecaca",
          borderRadius: 12,
          background: "#fef2f2",
          color: "#dc2626",
          fontSize: 13,
          padding: 16,
          textAlign: "center",
        }}
      >
        {errorMsg || "Failed to render molecule"}
      </div>
    );
  }

  return (
    <div style={{ border: "1px solid #e2e8f0", borderRadius: 12, padding: 8, background: "#fff", position: "relative" }}>
      {status === "loading" && (
        <div
          style={{
            position: "absolute",
            inset: 0,
            display: "flex",
            alignItems: "center",
            justifyContent: "center",
            background: "rgba(255,255,255,0.8)",
            borderRadius: 12,
            color: "#64748b",
            fontSize: 13,
          }}
        >
          Loading...
        </div>
      )}
      {svgContent && (
        <img
          src={svgContent}
          alt="Molecule structure"
          width={width}
          height={height}
          style={{ display: "block" }}
        />
      )}
    </div>
  );
}
