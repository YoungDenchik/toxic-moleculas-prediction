import { useEffect, useRef } from "react";

declare global {
  interface Window {
    SmilesDrawer: {
      Drawer: new (options: Record<string, unknown>) => {
        draw: (
          smiles: string,
          canvas: HTMLCanvasElement,
          theme: string,
          success?: () => void,
          error?: (err: unknown) => void
        ) => void;
      };
      parse: (smiles: string, callback: (tree: unknown) => void, errorCallback?: (err: unknown) => void) => void;
    };
  }
}

interface MoleculeViewerProps {
  smiles: string;
  width?: number;
  height?: number;
}

export function MoleculeViewer({ smiles, width = 300, height = 200 }: MoleculeViewerProps) {
  const canvasRef = useRef<HTMLCanvasElement>(null);

  useEffect(() => {
    if (!smiles || !canvasRef.current) return;

    const loadAndDraw = async () => {
      // Load SmilesDrawer if not already loaded
      if (!window.SmilesDrawer) {
        const script = document.createElement("script");
        script.src = "https://unpkg.com/smiles-drawer@2.1.7/dist/smiles-drawer.min.js";
        script.async = true;
        await new Promise<void>((resolve, reject) => {
          script.onload = () => resolve();
          script.onerror = () => reject(new Error("Failed to load SmilesDrawer"));
          document.head.appendChild(script);
        });
      }

      const drawer = new window.SmilesDrawer.Drawer({
        width,
        height,
        bondThickness: 1.5,
        bondLength: 20,
        shortBondLength: 0.85,
        bondSpacing: 0.18 * 20,
        atomVisualization: "default",
        isomeric: true,
        debug: false,
        terminalCarbons: false,
        explicitHydrogens: false,
        overlapSensitivity: 0.42,
        overlapResolutionIterations: 1,
        compactDrawing: false,
        fontSizeLarge: 11,
        fontSizeSmall: 6,
        padding: 20,
        experimentalSSSR: true,
        kkThreshold: 0.1,
        kkInnerThreshold: 0.1,
        kkMaxIteration: 20000,
        kkMaxInnerIteration: 50,
        kkMaxEnergy: 1e9,
      });

      window.SmilesDrawer.parse(
        smiles,
        (tree: unknown) => {
          drawer.draw(tree as unknown as string, canvasRef.current!, "light");
        },
        (err: unknown) => {
          console.error("Failed to parse SMILES:", err);
          const ctx = canvasRef.current?.getContext("2d");
          if (ctx) {
            ctx.clearRect(0, 0, width, height);
            ctx.fillStyle = "#ef4444";
            ctx.font = "14px system-ui";
            ctx.textAlign = "center";
            ctx.fillText("Invalid SMILES structure", width / 2, height / 2);
          }
        }
      );
    };

    loadAndDraw();
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

  return (
    <div style={{ border: "1px solid #e2e8f0", borderRadius: 12, padding: 8, background: "#fff" }}>
      <canvas ref={canvasRef} width={width} height={height} style={{ display: "block" }} />
    </div>
  );
}
