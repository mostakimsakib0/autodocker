module.exports = {
  hooks: {
    readPackage(pkg) {
      if (pkg.name === "ngl") {
        pkg.dependencies = {};
      }
      return pkg;
    }
  }
};
