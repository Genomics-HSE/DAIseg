{
  description = "daiseg dev environment";
  nixConfig.bash-prompt = "daiseg";

  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixpkgs-unstable";
    utils = {
      url = "github:numtide/flake-utils";
    };
  };
  outputs =
    {
      self,
      nixpkgs,
      utils,
    }:
    utils.lib.eachDefaultSystem (
      system:
      let
        pkgs = import nixpkgs { inherit system; };
        daiseg = pkgs.python3Packages.buildPythonApplication {
          pname = "daiseg";
          version = "0.1.0";
          pyproject = true;
          src = self;
          build-system = [ pkgs.python3Packages.setuptools ];

          dependencies = with pkgs.python3Packages; [
            numpy
            scipy
            pandas
            numba
            pysam
            jq
          ];

          # add stuff we call from subprocess to PATH
          makeWrapperArgs = [
            "--prefix PATH : ${
              pkgs.lib.makeBinPath (
                with pkgs;
                [
                  bcftools
                  bedtools
                  jq
                  gzip
                  findutils
                  gnused
                  coreutils
                ]
              )
            }"
          ];

          pythonImportsCheck = [ "daiseg.daiseg" ];
          doCheck = false;

          meta = {
            description = "DAIseg archaic segment inference";
            homepage = "https://github.com/Genomics-HSE/DAIseg";
            mainProgram = "daiseg";
            platforms = pkgs.lib.platforms.unix;
          };
        };
      in
      {

        packages = {
          inherit daiseg;
          default = daiseg;
        };

        apps.default = {
          type = "app";
          program = "${daiseg}/bin/daiseg";
        };

        devShells.default =
          with pkgs;
          mkShell {
            name = "daiseg";
            buildInputs = [
              wget
              ruff

              jq
              bedtools
              bcftools

              # htslib stuff for pysam
              zlib
              bzip2
              xz
              htslib
            ];

            shellHook = ''
              ${pkgs.kittysay}/bin/kittysay -t "mreow"
            '';
          };
      }
    );
}
