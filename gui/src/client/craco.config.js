module.exports = {
  webpack: {
    configure: (config) => {
      // Remove both the TypeScript checker and CRA’s ESLint plugin just in case
      const kill = new Set([
        'ForkTsCheckerWebpackPlugin', // TS type checker
        'ESLintPlugin',               // eslint-webpack-plugin
        'ESLintWebpackPlugin'         // alt name on some setups
      ]);

      // Log plugins so you can confirm removal the first time
      console.log('CRA plugins:', (config.plugins || []).map(p => p?.constructor?.name));

      config.plugins = (config.plugins || []).filter(
        (p) => !kill.has(p?.constructor?.name)
      );
      return config;
    },
  },
};