import { useEffect, useRef, useCallback, useState } from "react";
import { StandaloneStructServiceProvider } from "ketcher-standalone";
import { Editor } from "ketcher-react";
import type { Ketcher } from "ketcher-core";

import "ketcher-react/dist/index.css";

interface KetcherEditorProps {
  onSmilesChange: (smiles: string) => void;
  initialSmiles?: string;
}

const structServiceProvider = new StandaloneStructServiceProvider();

export function KetcherEditor({ onSmilesChange, initialSmiles }: KetcherEditorProps) {
  const ketcherRef = useRef<Ketcher | null>(null);
  const [isReady, setIsReady] = useState(false);
  const [error, setError] = useState<string | null>(null);
  const lastSmilesRef = useRef<string>("");
  const isSettingSmilesRef = useRef(false);

  const handleInit = useCallback((ketcher: Ketcher) => {
    console.log("[KetcherEditor] Initialized");
    ketcherRef.current = ketcher;
    setIsReady(true);
    setError(null);

    // Set initial SMILES if provided
    if (initialSmiles && initialSmiles.trim()) {
      isSettingSmilesRef.current = true;
      ketcher.setMolecule(initialSmiles).finally(() => {
        isSettingSmilesRef.current = false;
        lastSmilesRef.current = initialSmiles;
      });
    }
  }, [initialSmiles]);

  // Subscribe to structure changes
  useEffect(() => {
    if (!ketcherRef.current || !isReady) return;

    const ketcher = ketcherRef.current;

    const checkForChanges = async () => {
      if (isSettingSmilesRef.current) return;

      try {
        const smiles = await ketcher.getSmiles();
        if (smiles !== lastSmilesRef.current) {
          lastSmilesRef.current = smiles;
          onSmilesChange(smiles);
        }
      } catch {
        // Ignore errors when getting SMILES (e.g., empty structure)
      }
    };

    // Poll for changes since Ketcher doesn't have a reliable change event in standalone mode
    const pollInterval = setInterval(checkForChanges, 500);

    return () => {
      clearInterval(pollInterval);
    };
  }, [isReady, onSmilesChange]);

  // Handle external SMILES updates
  useEffect(() => {
    if (!ketcherRef.current || !isReady || !initialSmiles) return;
    if (initialSmiles === lastSmilesRef.current) return;

    isSettingSmilesRef.current = true;
    ketcherRef.current.setMolecule(initialSmiles).finally(() => {
      isSettingSmilesRef.current = false;
      lastSmilesRef.current = initialSmiles;
    });
  }, [initialSmiles, isReady]);

  const handleError = useCallback((message: string) => {
    console.error("[KetcherEditor] Error:", message);
    setError(message);
  }, []);

  return (
    <div style={{ position: "relative" }}>
      {!isReady && !error && (
        <div
          style={{
            position: "absolute",
            inset: 0,
            display: "flex",
            alignItems: "center",
            justifyContent: "center",
            background: "rgba(255,255,255,0.9)",
            zIndex: 10,
            borderRadius: 12,
          }}
        >
          Loading molecule editor...
        </div>
      )}
      {error && (
        <div
          style={{
            padding: 16,
            background: "#fef2f2",
            border: "1px solid #fecaca",
            borderRadius: 12,
            color: "#dc2626",
            marginBottom: 8,
          }}
        >
          Editor error: {error}
        </div>
      )}
      <div
        style={{
          height: 450,
          border: "1px solid #e2e8f0",
          borderRadius: 12,
          overflow: "hidden",
          background: "#fff",
        }}
      >
        <Editor
          staticResourcesUrl="."
          structServiceProvider={structServiceProvider}
          onInit={handleInit}
          errorHandler={handleError}
        />
      </div>
      <div style={{ marginTop: 8, fontSize: 12, color: "#64748b" }}>
        Draw a molecule structure. Changes are detected automatically.
      </div>
    </div>
  );
}
