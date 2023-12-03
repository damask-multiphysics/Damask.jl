using Documenter, Damask

format = Documenter.HTML(sidebar_sitename=false)
makedocs(sitename="Damask.jl",remotes=nothing,format=format)
