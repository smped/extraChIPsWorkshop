# extraChIPs Workshop

This package is taken from the template for building a Bioconductor workshop. The package
includes Github actions to:

1. Set up bioconductor/bioconductor_docker:devel on Github resources
2. Install package dependencies for your package (based on the `DESCRIPTION` file)
3. Run `rcmdcheck::rcmdcheck`
4. Build a pkgdown website and push it to github pages
5. Build a docker image with the installed package and dependencies and deploy to [the Github Container Repository](https://docs.github.com/en/packages/working-with-a-github-packages-registry/working-with-the-container-registry#pulling-container-images) at the name `ghcr.io/gihub_user/repo_name`, all lowercase. 


## To use the resulting image:

To try with **this** repository docker image, choose any password you'd like then:

```sh
docker run -e PASSWORD=<my_chosen_password> -p 8787:8787 ghcr.io/smped/extrachipsworkshop
```

*NOTE*: Running docker that uses the password in plain text like above exposes the password to others 
in a multi-user system (like a shared workstation or compute node). In practice, consider using an environment 
variable instead of plain text to pass along passwords and other secrets in docker command lines. 


