import { useState, useCallback, useEffect, useRef } from "react";
import { apiAnalyze } from "../api/client";
import type { AnalyzeResponse, UmapPoint, ChemicalSpace } from "../types";
import { MoleculeViewer } from "./MoleculeViewer";
import { LipinskiProperties } from "./LipinskiProperties";
import { ToxicityExplanation } from "./ToxicityExplanation";
import { UmapPlot } from "./UmapPlot";
import { KetcherEditor } from "./KetcherEditor";

interface MoleculeAnalyzerProps {
  backgroundPoints: UmapPoint[];
}

const EXAMPLE_MOLECULES = [
  { name: "Aspirin", smiles: "CC(=O)OC1=CC=CC=C1C(=O)O" },
  { name: "Caffeine", smiles: "CN1C=NC2=C1C(=O)N(C(=O)N2C)C" },
  { name: "Ibuprofen", smiles: "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O" },
  { name: "Ethanol", smiles: "CCO" },
  { name: "Dofetilide", smiles: "CC(C)CCNC(=O)c1cc(ccc1OC)NS(=O)(=O)c2ccc(cc2)NCCCC(C)C" },
  { name: "Terfenadine", smiles: "CC(C)c1ccc(cc1)C(O)CCCN2CCC(CC2)C(O)(c3ccccc3)c4ccccc4" },
];

export function MoleculeAnalyzer({ backgroundPoints }: MoleculeAnalyzerProps) {
  const [smiles, setSmiles] = useState("");
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState<string | null>(null);
  const [result, setResult] = useState<AnalyzeResponse | null>(null);
  const [inputMode, setInputMode] = useState<"text" | "draw">("text");
  const debounceRef = useRef<ReturnType<typeof setTimeout> | null>(null);

  const analyzeSmiles = useCallback(async (smilesInput: string) => {
    const trimmed = smilesInput.trim();
    if (!trimmed) {
      setResult(null);
      setError(null);
      return;
    }

    setLoading(true);
    setError(null);

    try {
      const response = await apiAnalyze(trimmed);
      setResult(response);
      setError(null);
    } catch (err) {
      setError(err instanceof Error ? err.message : "Failed to analyze molecule");
      setResult(null);
    } finally {
      setLoading(false);
    }
  }, []);

  // Debounced analysis for drawing mode
  const debouncedAnalyze = useCallback(
    (smilesInput: string) => {
      if (debounceRef.current) {
        clearTimeout(debounceRef.current);
      }
      debounceRef.current = setTimeout(() => {
        analyzeSmiles(smilesInput);
      }, 800); // Slightly longer debounce for drawing mode
    },
    [analyzeSmiles]
  );

  // Handle drawing mode SMILES changes
  const handleDrawingChange = useCallback(
    (newSmiles: string) => {
      setSmiles(newSmiles);
      if (newSmiles.trim()) {
        debouncedAnalyze(newSmiles);
      } else {
        setResult(null);
        setError(null);
      }
    },
    [debouncedAnalyze]
  );

  // Cleanup debounce on unmount
  useEffect(() => {
    return () => {
      if (debounceRef.current) {
        clearTimeout(debounceRef.current);
      }
    };
  }, []);

  const handlePredict = () => {
    analyzeSmiles(smiles);
  };

  const handleClear = () => {
    setSmiles("");
    setResult(null);
    setError(null);
  };

  const handleExampleClick = (exampleSmiles: string) => {
    setSmiles(exampleSmiles);
    analyzeSmiles(exampleSmiles);
  };

  const userMoleculePosition: { position: ChemicalSpace; smiles: string; isToxic: boolean } | undefined = result
    ? {
        position: result.chemical_space,
        smiles: result.canonical_smiles,
        isToxic: result.prediction.label === "toxic",
      }
    : undefined;

  return (
    <div style={{ display: "flex", flexDirection: "column", gap: 24 }}>
      {/* Mode Toggle */}
      <div style={{ display: "flex", gap: 8 }}>
        <button
          onClick={() => setInputMode("text")}
          style={{
            padding: "10px 20px",
            borderRadius: 10,
            border: "2px solid",
            borderColor: inputMode === "text" ? "#1d4ed8" : "#e2e8f0",
            background: inputMode === "text" ? "#1d4ed8" : "#fff",
            color: inputMode === "text" ? "#fff" : "#0f172a",
            fontWeight: 600,
            cursor: "pointer",
            transition: "all 0.2s",
            display: "flex",
            alignItems: "center",
            gap: 8,
          }}
        >
          <span style={{ fontSize: 18 }}>📝</span>
          SMILES Input
        </button>
        <button
          onClick={() => setInputMode("draw")}
          style={{
            padding: "10px 20px",
            borderRadius: 10,
            border: "2px solid",
            borderColor: inputMode === "draw" ? "#1d4ed8" : "#e2e8f0",
            background: inputMode === "draw" ? "#1d4ed8" : "#fff",
            color: inputMode === "draw" ? "#fff" : "#0f172a",
            fontWeight: 600,
            cursor: "pointer",
            transition: "all 0.2s",
            display: "flex",
            alignItems: "center",
            gap: 8,
          }}
        >
          <span style={{ fontSize: 18 }}>🎨</span>
          Draw Molecule
        </button>
      </div>

      {/* SMILES Text Input Mode */}
      {inputMode === "text" && (
        <div style={{ display: "grid", gridTemplateColumns: "1fr 1fr", gap: 20 }}>
          {/* Left: Input */}
          <div style={{ display: "flex", flexDirection: "column", gap: 16 }}>
            <div style={{ border: "1px solid #e2e8f0", borderRadius: 16, padding: 16 }}>
              <label style={{ display: "block", fontSize: 13, color: "#64748b", marginBottom: 8, fontWeight: 600 }}>
                Enter SMILES notation
              </label>
              <textarea
                value={smiles}
                onChange={(e) => setSmiles(e.target.value)}
                placeholder='e.g., "CCO" for ethanol or "c1ccccc1" for benzene'
                rows={4}
                style={{
                  width: "100%",
                  boxSizing: "border-box",
                  resize: "vertical",
                  border: "1px solid #cbd5e1",
                  borderRadius: 12,
                  padding: 12,
                  fontFamily: "ui-monospace, SFMono-Regular, Menlo, Monaco, Consolas, monospace",
                  fontSize: 14,
                  outline: "none",
                }}
              />

              <div style={{ display: "flex", gap: 10, marginTop: 12, alignItems: "center" }}>
                <button
                  onClick={handlePredict}
                  disabled={loading || !smiles.trim()}
                  style={{
                    background: loading || !smiles.trim() ? "#94a3b8" : "#1d4ed8",
                    color: "white",
                    padding: "12px 24px",
                    borderRadius: 12,
                    border: "none",
                    fontWeight: 700,
                    fontSize: 14,
                    cursor: loading || !smiles.trim() ? "not-allowed" : "pointer",
                    display: "flex",
                    alignItems: "center",
                    gap: 8,
                  }}
                >
                  {loading && <Spinner size={14} />}
                  {loading ? "Analyzing..." : "Analyze Molecule"}
                </button>

                <button
                  onClick={handleClear}
                  style={{
                    background: "white",
                    color: "#0f172a",
                    padding: "12px 24px",
                    borderRadius: 12,
                    border: "1px solid #cbd5e1",
                    fontWeight: 700,
                    fontSize: 14,
                    cursor: "pointer",
                  }}
                >
                  Clear
                </button>
              </div>

              {/* Example molecules */}
              <div style={{ marginTop: 16 }}>
                <div style={{ fontSize: 12, color: "#64748b", marginBottom: 8 }}>Try an example:</div>
                <div style={{ display: "flex", flexWrap: "wrap", gap: 6 }}>
                  {EXAMPLE_MOLECULES.map((mol) => (
                    <button
                      key={mol.name}
                      onClick={() => handleExampleClick(mol.smiles)}
                      style={{
                        padding: "6px 12px",
                        borderRadius: 999,
                        border: "1px solid #e2e8f0",
                        background: "#f8fafc",
                        fontSize: 12,
                        fontWeight: 500,
                        cursor: "pointer",
                        color: "#334155",
                        transition: "all 0.15s",
                      }}
                      onMouseOver={(e) => {
                        e.currentTarget.style.background = "#e2e8f0";
                      }}
                      onMouseOut={(e) => {
                        e.currentTarget.style.background = "#f8fafc";
                      }}
                    >
                      {mol.name}
                    </button>
                  ))}
                </div>
              </div>
            </div>

            {/* Error display */}
            {error && (
              <div
                style={{
                  padding: "12px 16px",
                  background: "#fef2f2",
                  border: "1px solid #fecaca",
                  borderRadius: 12,
                  color: "#991b1b",
                  fontSize: 13,
                }}
              >
                <strong>Error:</strong> {error}
              </div>
            )}
          </div>

          {/* Right: Molecule Visualization */}
          <div style={{ display: "flex", flexDirection: "column", gap: 16 }}>
            <div style={{ border: "1px solid #e2e8f0", borderRadius: 16, padding: 16 }}>
              <div style={{ fontWeight: 800, fontSize: 14, marginBottom: 12 }}>Molecule Structure</div>
              <MoleculeViewer smiles={result?.canonical_smiles || smiles} width={350} height={250} />
              {result?.canonical_smiles && (
                <div style={{ marginTop: 8, fontSize: 11, color: "#64748b" }}>
                  Canonical: <code style={{ background: "#f1f5f9", padding: "2px 6px", borderRadius: 4 }}>{result.canonical_smiles}</code>
                </div>
              )}
            </div>
          </div>
        </div>
      )}

      {/* Drawing Mode with Ketcher */}
      {inputMode === "draw" && (
        <div style={{ display: "flex", flexDirection: "column", gap: 16 }}>
          <div style={{ display: "grid", gridTemplateColumns: "2fr 1fr", gap: 20 }}>
            {/* Ketcher Editor */}
            <div style={{ border: "1px solid #e2e8f0", borderRadius: 16, padding: 16 }}>
              <div style={{ fontWeight: 800, fontSize: 14, marginBottom: 12 }}>
                Draw Your Molecule
              </div>
              <KetcherEditor onSmilesChange={handleDrawingChange} initialSmiles={smiles} />
              {smiles && (
                <div style={{ marginTop: 12, padding: 10, background: "#f1f5f9", borderRadius: 8 }}>
                  <div style={{ fontSize: 11, color: "#64748b", marginBottom: 4 }}>Generated SMILES:</div>
                  <code style={{ fontSize: 12, color: "#0f172a", wordBreak: "break-all" }}>{smiles}</code>
                </div>
              )}
              {loading && (
                <div style={{ marginTop: 12, display: "flex", alignItems: "center", gap: 8, color: "#1d4ed8" }}>
                  <Spinner size={16} />
                  <span style={{ fontSize: 13, fontWeight: 500 }}>Analyzing...</span>
                </div>
              )}
            </div>

            {/* Molecule Visualization */}
            <div style={{ display: "flex", flexDirection: "column", gap: 16 }}>
              <div style={{ border: "1px solid #e2e8f0", borderRadius: 16, padding: 16 }}>
                <div style={{ fontWeight: 800, fontSize: 14, marginBottom: 12 }}>Structure Preview</div>
                <MoleculeViewer smiles={result?.canonical_smiles || smiles} width={280} height={200} />
              </div>

              {/* Error display */}
              {error && (
                <div
                  style={{
                    padding: "12px 16px",
                    background: "#fef2f2",
                    border: "1px solid #fecaca",
                    borderRadius: 12,
                    color: "#991b1b",
                    fontSize: 13,
                  }}
                >
                  <strong>Error:</strong> {error}
                </div>
              )}
            </div>
          </div>
        </div>
      )}

      {/* Results Section */}
      {result && (
        <div style={{ display: "grid", gridTemplateColumns: "1fr 1fr", gap: 20 }}>
          {/* Left column: Prediction and Lipinski */}
          <div style={{ display: "flex", flexDirection: "column", gap: 16 }}>
            <ToxicityExplanation prediction={result.prediction} explainability={result.explainability} />
            <LipinskiProperties properties={result.properties} lipinskiPass={result.lipinski_pass} />
          </div>

          {/* Right column: UMAP */}
          <UmapPlot points={backgroundPoints} userMolecule={userMoleculePosition} />
        </div>
      )}

      {/* Loading overlay for results */}
      {loading && result && (
        <div
          style={{
            position: "fixed",
            top: 0,
            left: 0,
            right: 0,
            bottom: 0,
            background: "rgba(255,255,255,0.7)",
            display: "flex",
            alignItems: "center",
            justifyContent: "center",
            zIndex: 100,
          }}
        >
          <div
            style={{
              padding: "20px 40px",
              background: "#fff",
              borderRadius: 16,
              boxShadow: "0 4px 20px rgba(0,0,0,0.1)",
              display: "flex",
              alignItems: "center",
              gap: 12,
            }}
          >
            <Spinner size={24} />
            <span style={{ fontWeight: 600 }}>Analyzing molecule...</span>
          </div>
        </div>
      )}

      {/* Spinner animation */}
      <style>{`
        @keyframes spin {
          to { transform: rotate(360deg); }
        }
      `}</style>
    </div>
  );
}

function Spinner({ size = 16 }: { size?: number }) {
  return (
    <span
      style={{
        width: size,
        height: size,
        border: `${Math.max(2, size / 8)}px solid #1d4ed8`,
        borderTopColor: "transparent",
        borderRadius: "50%",
        animation: "spin 1s linear infinite",
        display: "inline-block",
      }}
    />
  );
}
