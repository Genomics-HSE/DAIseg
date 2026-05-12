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
      in
      {
        devShells.default =
          with pkgs;
          mkShell {
            name = "daiseg";
            buildInputs = [
              wget
              ouch
              ruff
              jq
              bedtools
              bcftools
            ];

            shellHook = ''
              ${pkgs.kittysay}/bin/kittysay -t "mreow"
            '';
          };
      }
    );
}
