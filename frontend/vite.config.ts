import { defineConfig } from 'vite'
import react from '@vitejs/plugin-react'
import { nodePolyfills } from 'vite-plugin-node-polyfills'

// https://vite.dev/config/
export default defineConfig({
  plugins: [
    react(),
    nodePolyfills({
      // Enable polyfills for specific globals
      globals: {
        Buffer: true,
        global: true,
        process: true,
      },
      // Enable polyfills for specific modules
      protocolImports: true,
    }),
  ],
  define: {
    // Fix for some packages that check for process.env
    'process.env': {},
  },
  optimizeDeps: {
    // Include ketcher packages for optimization
    include: ['ketcher-react', 'ketcher-core', 'ketcher-standalone'],
  },
})
