using Pkg
if isfile("Project.toml") && isfile("Manifest.toml")
    Pkg.activate(".")

    Pkg.instantiate()
end
using Revise
import PeriLab
using TimerOutputs
using LoggingExtras
const to = TimerOutput()

dry_run = false
verbose = false
debug = false
filename = ""
for i in eachindex(ARGS)
    arg = ARGS[i]

    if arg == "--dry_run"
        global dry_run = true
    elseif arg == "-v" || arg == "--verbose"
        global verbose = true
    elseif arg == "-d" || arg == "--debug"
        global debug = true
    else
        global filename = arg
    end
end

if debug
    demux_logger = TeeLogger(
        MinLevelLogger(FileLogger(split(filename, ".")[1] * ".log"), Logging.Debug),
        MinLevelLogger(ConsoleLogger(stderr), Logging.Debug),
    )
else
    demux_logger = TeeLogger(
        MinLevelLogger(FileLogger(split(filename, ".")[1] * ".log"), Logging.Info),
        MinLevelLogger(ConsoleLogger(stderr), Logging.Info),
    )
end
global_logger(demux_logger)

@timeit to "PeriLab" begin
    PeriLab.main(filename, to, dry_run, verbose)
end

if verbose
    @info demux_logger to
end